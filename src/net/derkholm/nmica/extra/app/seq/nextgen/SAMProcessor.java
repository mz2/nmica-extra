package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.DummySymbolList;
import org.biojava.bio.symbol.RangeLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@App(overview = "Create an SQLite database usable by nextgenseq-dazzle", generateStub = true)
@NMExtraApp(launchName = "ngsamprocess", vm = VirtualMachine.SERVER)
public class SAMProcessor {
	private SAMFileReader inReader;
	private SequenceDB seqDB = new HashSequenceDB();
	private int qualityCutoff = 10;
	private int expandReadsBy;
	private File indexFile;
	private String in;
	private Map<String,Integer> refSeqLengths = new HashMap<String,Integer>();

	@Option(help="Input map in SAM/BAM format")
	public void setMap(String str) {
		System.err.println("Reading map...");
		this.in = str;
	}
	
	@Option(help="Index file for the reads")
	public void setIndex(File f) {
		this.indexFile = f;
	}
	
	
	@Option(help="Expand reads by specified number of nucleotides (bound by reference sequence ends)", optional=true)
	public void setExpandReadsBy(int i) {
		this.expandReadsBy = i;
	}

	
	@Option(help="Reference sequence names and lengths in a TSV formatted file")
	public void setRefLengths(File f) throws NoSuchElementException, BioException, NumberFormatException, IOException {
		BufferedReader reader = new BufferedReader(new FileReader(f));
		
		String line = null;
		while ((line = reader.readLine()) != null) {
			StringTokenizer tok = new StringTokenizer(line,"\t");
			refSeqLengths.put(tok.nextToken(), Integer.parseInt(tok.nextToken()));
		}
	}
	
	@Option(help="Mapping quality threshold (exclude reads whose mapping quality is below. default=10)", optional=true)
	public void setMappingQualityCutoff(int quality) {
		this.qualityCutoff = quality;
	}

	public void main(String[] args) throws BioException {
		this.inReader = new SAMFileReader(new File(in),indexFile);
		this.inReader.setValidationStringency(ValidationStringency.SILENT);
		
		int excludedReads = 0;
		int readCount = 0;
		Sequence seq = null;
		String source = "samprocessor";
		String type = "read";

		Set<String> seenRefSeqNames = new TreeSet<String>();

		//this.inReader.queryOverlapping(sequence, start, end);
		for (final SAMRecord samRecord : this.inReader) {
			readCount += 1;

			int quality = samRecord.getMappingQuality();
			if (quality < qualityCutoff) {
				excludedReads += 1;
				continue;
			}
			String refName = samRecord.getReferenceName();
			assert refName != null;
			seenRefSeqNames.add(refName);

			int start = samRecord.getAlignmentStart();
			int end = samRecord.getAlignmentEnd();


			boolean isPosStrand = !samRecord.getReadNegativeStrandFlag();

			if (isPosStrand) {
				end = Math.min(seq.length()-1, end + this.expandReadsBy);
			} else {
				start = Math.max(0, start - this.expandReadsBy);
			}

			if (seq == null) {
				System.err.println("Collecting reads mapped to reference sequence " + refName);
				
				DummySymbolList symList = new DummySymbolList(DNATools.getDNA(), Integer.MAX_VALUE);
				seq = new SimpleSequence(symList,refName,refName,Annotation.EMPTY_ANNOTATION);
			} else if (seq.getName().equals(refName)) {
				StrandedFeature.Template template = new StrandedFeature.Template();
				template.source = source;
				template.type = type;
				template.location = new RangeLocation(start, end);
				template.strand = isPosStrand ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE;

				seq.createFeature(template);
			} else {
				assert !seenRefSeqNames.contains(refName);

				this.processSequence(seq);
				
				DummySymbolList symList = new DummySymbolList(DNATools.getDNA(), Integer.MAX_VALUE);
				seq = new SimpleSequence(symList,refName,refName,Annotation.EMPTY_ANNOTATION);
				
				System.err.println("Collecting reads mapped to reference sequence " + refName);
			}
		}

		System.err.printf("Excluded %d reads (%.2f%)%n", excludedReads, (double)excludedReads / (double)readCount * 100.0);
	}

	private void processSequence(Sequence seq) {

		System.err.println("Processing "+ seq.getName());
	}
}
