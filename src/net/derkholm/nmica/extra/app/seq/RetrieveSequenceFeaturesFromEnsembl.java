package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.sql.SQLException;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Get noncoding sequences from Ensembl for motif discovery", generateStub = true)
@NMExtraApp(launchName = "nmensemblfeat", vm = VirtualMachine.SERVER)
public class RetrieveSequenceFeaturesFromEnsembl extends RetrieveEnsemblSequences {

	private File outFile;
	private File featuresFile;
	private int expandToLength = 0;
	private int minNonN = 1;
	private int maxDistFromGene;

	@Option(help="Output file",optional=true)
	public void setOut(File f) {
		this.outFile = f;
	}

	@Option(help="Features file (read from stdin if not included)", optional=true)
	public void setFeatures(File f) {
		this.featuresFile = f;
	}
	
	@Option(help="Expand to length", optional=true)
	public void setExpandToLength(int i) {
		this.expandToLength = i;
	}
	
	@Option(help="Minimun number of gap symbols (N)", optional=true)
	public void setMinNonN(int i) {
		this.minNonN = i;
	}
	
	@Option(help="Maximum distance from nearby gene (do not output if not found)", optional=true)
	public void setMaxDistanceFromGene(int i) {
		this.maxDistFromGene = i;
	}
	

	public static int gapSymbolCount(SymbolList seq) {
		int numNs = 0;
		for (Iterator<?> i = seq.iterator(); i.hasNext(); ) {
            Symbol s = (Symbol) i.next();
            if (s == DNATools.n() || s == seq.getAlphabet().getGapSymbol()) {
                ++numNs;
            }
        }
		
		return numNs;
	}

	
	public void main(String[] args) throws SQLException, Exception {
		initializeEnsemblConnection();

		final OutputStream os;
		if (this.outFile == null) {
			os = System.out;
		} else {
			os = new FileOutputStream(this.outFile);
		}

		InputStream inputStream;
		if (featuresFile == null) {
			inputStream = System.in;
		} else {
			inputStream = new FileInputStream(this.featuresFile);
		}

		GFFParser parser = new GFFParser();
		parser.parse(
				new BufferedReader(new InputStreamReader(inputStream)),
				new GFFDocumentHandler() {
			public void commentLine(String str) {}
			public void endDocument() {}
			public void startDocument(String str) {}

			public void recordLine(GFFRecord recLine) {
				System.err.printf(".");
				try {

					int start = recLine.getStart();
					int end = recLine.getEnd();
					
					if (expandToLength > 0) {
						if ((end - start + 1) < expandToLength) {
							start = Math.max(1, start - (expandToLength / 2));
							end = end + (expandToLength / 2);
						}
					}
					
					if (maxDistFromGene > 0) {
						if (recLine.getStrand().equals(StrandedFeature.UNKNOWN)) {
							System.err.println("WARNING: Feature does not have strand. Cannot determine distance from genes on the same strand");
							return;
						}
						
						FeatureHolder transcripts = seqDB.filter(
														new FeatureFilter.And(
															new FeatureFilter.OverlapsLocation(
																new RangeLocation(recLine.getStart(),recLine.getEnd())),
															new FeatureFilter.BySequenceName(recLine.getSeqName())));
						
						transcripts = transcripts.filter(new FeatureFilter.StrandFilter(recLine.getStrand()));
						
						StrandedFeature closestTranscript = null;
						Set<StrandedFeature> features;
						
						if (recLine.getStrand().equals(StrandedFeature.POSITIVE)) {
							features = new TreeSet<StrandedFeature>(new DistanceComparator(recLine.getStart()));
						} else if (recLine.getStrand().equals(StrandedFeature.NEGATIVE)) {
							features = new TreeSet<StrandedFeature>(new DistanceComparator(recLine.getEnd()));
						} else {
							throw new BioError("Feature must be stranded");
						}
						
						for (Iterator<?> fi = transcripts.features(); fi.hasNext();) {
							StrandedFeature transcript = (StrandedFeature) fi.next();
							Location loc = transcript.getLocation();
							int tStart,tEnd;
							
							if (recLine.getStrand().equals(StrandedFeature.POSITIVE)) {
								tStart = loc.getMin();
								tEnd = loc.getMax();
							} else {
								tStart = loc.getMax();
								tEnd = loc.getMin();
							}

							transcript.getAnnotation().getProperty("ensembl.xrefs");
							
							Object o = transcript.getAnnotation().getProperty("ensembl.xrefs");
							if (ignoreGenesWithNoCrossReferences && 
									((o == null) || 
									(!(o instanceof List)) || 
									(((List)o)).size() == 0)) {
								System.err.printf("Transcript of gene %s does not have cross-references%n", transcript.getAnnotation().getProperty("ensembl.id"));
							}
							
							RangeLocation tssLocationRange = new RangeLocation(tStart, tEnd);
							RangeLocation posLocationRange = new RangeLocation(
									recLine.getStart()-maxDistFromGene,
									recLine.getEnd()+maxDistFromGene);
							
							if (tssLocationRange.overlaps(posLocationRange)) {
								System.err.println(String.format(
									"Feature is within +/- %d the transcript of gene %s",
									maxDistFromGene,
									transcript.getAnnotation().getProperty("ensembl.gene_id")));
								features.add(transcript);
							}
						}
					}
					
					SymbolList symList = 
						seqDB.getSequence(
							recLine.getSeqName()).subList(start,end);
					
					if (recLine.getStrand().equals(StrandedFeature.NEGATIVE)) {
						symList = DNATools.reverseComplement(symList);
					}
					
					if (minNonN > 0 && 
						RetrieveSequenceFeaturesFromEnsembl
							.gapSymbolCount(symList) > minNonN) {
						return;
					}
					
					Sequence s = new SimpleSequence(symList, null,
							String.format("%s;%d-%d(%s)",
									recLine.getSeqName(),
									Math.max(1,start),
									end,
									recLine.getStrand()
										.equals(StrandedFeature.POSITIVE)? "+" : "-"),
							Annotation.EMPTY_ANNOTATION);

					RichSequence.IOTools.writeFasta(os, s, null);

				} catch (IllegalIDException e) {
					e.printStackTrace();
				} catch (BioException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		});
		System.err.println();
	}
	
	private class DistanceComparator implements Comparator<StrandedFeature> {
		private int referencePosition;
		
		public DistanceComparator(int refPos) {
			this.referencePosition = refPos;
		}
		
		public int compare(StrandedFeature o1, StrandedFeature o2) {
			int ref0Start, ref0End;
			ref0Start = Math.abs(o1.getLocation().getMin() - referencePosition);
			ref0End = Math.abs(o2.getLocation().getMax() - referencePosition);
			int ref0 = Math.min(ref0Start, ref0End);
			
			int ref1Start, ref1End;
			ref1Start = Math.abs(o1.getLocation().getMin() - referencePosition);
			ref1End = Math.abs(o1.getLocation().getMax() - referencePosition);
			int ref1 = Math.min(ref1Start,ref1End);
			
			if (ref0 < ref1) {
				return 1;
			} else if (ref0 > ref1) {
				return -1;
			} 
			
			return 0;			
		}
	}
}