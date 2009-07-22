package net.derkholm.nmica.extra.app;

import gfftools.GFFUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Write out sequence regions covered by GFF features", generateStub = true)
@NMExtraApp(launchName = "nmcoveredseq", vm = VirtualMachine.SERVER)
public class WriteCoveredSequences {
    public static enum Format {
    	GFF,
    	FASTA,
    	TABLE
    }
	
	private boolean negate;
	private File gffFile;
	private File seqFile;
	private Format outputFormat = Format.FASTA;

	@Option(help="Input feature file")
	public void setFeatures(File f) {
		this.gffFile = f;
	}
	
	@Option(help="Input sequence file")
	public void setSeqs(File f) {
		this.seqFile = f;
	}
	
	@Option(help="Output format: fasta|gff", optional=true)
	public void setFormat(Format format) {
		this.outputFormat = format;
	}
	
	@Option(help="Negate the output, i.e. output sequences that are NOT covered by the features.", optional=true)
	public void setNegate(boolean b) {
		this.negate = b;
	}

	public void main(String[] args)
	throws Exception
	{   

		Map<String,Location> locs = GFFUtils.gffToLocationMap(gffFile);
		for (SequenceIterator si = RichSequence
									.IOTools
										.readFastaDNA(
											new BufferedReader(
												new FileReader(seqFile)), null); si.hasNext();) {
			Sequence seq = si.nextSequence();
			Location loc = locs.get(seq.getName());
			
			org.biojava.bio.seq.StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;

			
			GFFWriter writer = null;
			if (outputFormat == Format.GFF) {
				writer = new GFFWriter(new PrintWriter(System.out));
			}
			
			if (loc == null && negate) {
				if (outputFormat == Format.FASTA) {
					RichSequence.IOTools.writeFasta(System.out, seq, null);					
				} else {
					for (Iterator<?> li = loc.blockIterator(); li.hasNext(); ) {
						Location l = (Location)li.next();
						GFFRecord r = new SimpleGFFRecord(seq.getName(),"nmcoveredseq","covered",l.getMin(),l.getMax(),Double.NaN,strand,0,"",null);
						
						writer.recordLine(r);
					}
				}
			} else if (loc == null &! negate) {
				//do nothing (nothing was covered by a location)
			} else {
				Location wanted;
				if (negate) {
					wanted = LocationTools.subtract(new RangeLocation(1, seq.length()), loc);					
				} else {
					wanted = loc;
				}
				
				for (Iterator<?> bi = wanted.blockIterator(); bi.hasNext(); ) { 
					Location wl = (Location) bi.next();
					
					if (outputFormat == Format.FASTA) {
						RichSequence.IOTools.writeFasta(System.out, new SimpleSequence(
								seq.subList(wl.getMin(), wl.getMax()),
								null,
								String.format("%s_|%d-%d|", seq.getName(), wl.getMin(), wl.getMax()),
								Annotation.EMPTY_ANNOTATION
						), null);						
					} else {
						System.err.printf("%s %s %s %s%n",seq,seq.getName(),wl,strand);
						GFFRecord r = new SimpleGFFRecord(seq.getName(),"nmcoveredseq","covered",wl.getMin(),wl.getMax(),Double.NaN,strand,0,"",new TreeMap());
						writer.recordLine(r);
					}
				}
			}   
			
			if (writer != null) writer.endDocument();
		}   
	}   
}