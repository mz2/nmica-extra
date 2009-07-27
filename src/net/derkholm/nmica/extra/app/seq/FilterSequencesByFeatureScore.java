package net.derkholm.nmica.extra.app.seq;

import gfftools.GFFUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.NoSuchElementException;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

import biobits.utils.IOTools;

@App(overview = "Write out sequence regions covered by GFF features", generateStub = true)
@NMExtraApp(launchName = "nmfilterseq", vm = VirtualMachine.SERVER)
public class FilterSequencesByFeatureScore {

    public static enum Format {
    	GFF,
    	FASTA
    }
	
	private Format format = Format.GFF;
	
	private double minScore = Double.NaN;
	@Option(help="Minimum score", optional=true)
	public void setMinScore(double s) {
		this.minScore = s;
	}

	private double maxScore = Double.NaN;
	private File featureFile;
	private File seqFile;
	private File seqsFile;
	private HashSequenceDB seqDB;
	private boolean validate;

	private boolean negate;

	@Option(help="Maximum score", optional=true)
	public void setMaxScore(double s) {
		this.maxScore = s;
	}

	@Option(help="Validate input (check that sequence identifiers match)", optional=true, userLevel = UserLevel.DEBUG)
	public void setValidate(boolean b) {
		this.validate = b;
	}
	
	@Option(help="Negate output (output features/sequences that don't fill the filter constraints", optional=true)
	public void setNegate(boolean b) {
		this.negate = b;
	}

	@Option(help = "Input sequences (output the intersecting sequences rather than a new GFF file)", optional=true)
	public void setSeqs(File f) throws FileNotFoundException, ChangeVetoException, NoSuchElementException, BioException {
		this.seqsFile = f;

		this.seqDB = new HashSequenceDB();
		for (SequenceIterator si = RichSequence
				.IOTools
				.readFastaDNA(
						new BufferedReader(
								new FileReader(f)), null); si.hasNext();) {

			seqDB.addSequence(si.nextSequence());
		}
		this.format = Format.FASTA;
	}

	@Option(help="Input file for features")
	public void setFeatures(File file) {
		this.featureFile = file;
	}

	public void main(String[] args) throws Exception {

		if (validate && (seqDB != null)) {
			Map<String,Location> locs = GFFUtils.gffToLocationMap(featureFile);
			WriteCoveredSequences.validateGFFSequenceIdentifiersAgainstSequences(locs, seqDB);
		}
		
		GFFParser gffp = new GFFParser();
		gffp.parse(IOTools.fileBufferedReader(featureFile), 
				new GFFFilteringDocumentHandler(this.minScore,this.maxScore,this.negate,this.format, seqDB));
		
		
		
	}

	private static class GFFFilteringDocumentHandler implements GFFDocumentHandler {
		private double minScore = Double.NaN;
		private double maxScore = Double.NaN;
		private Format format;
		private GFFWriter gffWriter;
		private SequenceDB sequences;
		private boolean negate;

		public GFFFilteringDocumentHandler() {
			
		}
		
		private GFFFilteringDocumentHandler(double minScore, double maxScore, boolean negate, Format format, SequenceDB sequences) {
			this.minScore = minScore;
			this.maxScore = maxScore;
			this.negate = negate;
			this.gffWriter = new GFFWriter(new PrintWriter(System.out));
			this.format = format;
			this.sequences = sequences;
		}
		
		public void startDocument(String locator) {
		}

		public void endDocument() {
		}

		public void commentLine(String comment) {
		}

		public void recordLine(GFFRecord record) {
			boolean allowOutput = true;
			if (!Double.isNaN(this.minScore)) {
				if (this.minScore > record.getScore()) {
					allowOutput = false;
				}
			}
			
			if (!Double.isNaN(this.maxScore)) {
				if (this.maxScore < record.getScore()) {
					allowOutput = false;
				}
			}
			
			if (negate) {
				allowOutput = !allowOutput;
			}
			
			if (allowOutput) {
				if (this.format == Format.GFF) {
					this.gffWriter.recordLine(record);					
				} else {
					Sequence seq;
					try {
						seq = new SimpleSequence(sequences.getSequence(record.getSeqName()).subList(record.getStart(),record.getEnd()),null,
								String.format("%s_|%d-%d|",record.getSeqName(),record.getStart(),record.getEnd()),Annotation.EMPTY_ANNOTATION);
					} catch (IllegalIDException e) {
						throw new BioError(e);
					} catch (BioException e) {
						throw new BioError(e);
					}
                	try {
						RichSequence.IOTools.writeFasta(System.out, seq, null);
					} catch (IOException e) {
						throw new BioError(e);
					}
				}
			}
		}
	}
}