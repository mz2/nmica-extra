package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedSet;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakEntry;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakFormat;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.RankOrder;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Convert peak files from various formats to GFF", generateStub = true)
@NMExtraApp(launchName = "ngpeak2gff", vm = VirtualMachine.SERVER)
public class PeaksToGFF {
	
	private PeakFormat format;
	private FileReader peaksReader;
	private RankOrder rankOrder = RetrievePeakSequencesFromEnsembl.RankOrder.DESC;
	private int aroundPeak;
	private int maxLength;
	private int minLength;
	private int maxCount;

	@Option(help="Peak file format")
	public void setFormat(RetrievePeakSequencesFromEnsembl.PeakFormat format) {
		this.format = format;
	}

	@Option(help="Input peak file")
	public void setPeaks(FileReader f) {
		this.peaksReader = f;
	}
	
	@Option(help="Positions to fetch around peak maximae", optional=true)
	public void setAroundPeak(int i) {
		this.aroundPeak = i;
	}
	
	@Option(help="Maximum peak length (not defined by default)",optional=true)
	public void setMaxLength(int maxLength) {
		this.maxLength = maxLength;
	}
	
	@Option(help="Minimum peak length (default=20)",optional=true)
	public void setMinLength(int minLength) {
		this.minLength = minLength;
	}
	
	@Option(help="Maximum number of peaks to output", optional=true)
	public void setMaxCount(int maxCount) {
		this.maxCount = maxCount;
	}
	
	
	public void main(String[] args) throws FileNotFoundException, IOException {
		SortedSet<PeakEntry> peaks = RetrievePeakSequencesFromEnsembl.parsePeaks(
				new BufferedReader(peaksReader), 
				format, 
				rankOrder, 
				aroundPeak, 
				minLength, 
				maxLength);
		
		if (maxCount == 0) {
			maxCount = peaks.size();
		}
		
		System.err.printf("Parsed %d peaks", peaks.size());
		GFFWriter writer = new GFFWriter(new PrintWriter(System.out));
		
		Iterator<PeakEntry> peakIterator = peaks.iterator();
		
		int i = 0;
		while (i < maxCount) {
			final PeakEntry peak = peakIterator.next();
			
			writer.recordLine(new GFFRecord() {

				public String getComment() {
					return String.format("\tp-value:%.3f fdr:%.3f", peak.pValue,peak.fdr);
				}

				public int getEnd() {
					return peak.endCoord;
				}

				public String getFeature() {
					return "peak";
				}

				public int getFrame() {
					return 0;
				}

				public Map<String, Object> getGroupAttributes() {
					Map<String, Object> map = new HashMap<String,Object>();
					//map.put("p-value", ""+peak.pValue);
					return map;
				}

				public double getScore() {
					return peak.foldChange;
				}

				public String getSeqName() {
					return peak.seqName;
				}

				public String getSource() {
					return format.name();
				}

				public int getStart() {
					return peak.startCoord;
				}

				public Strand getStrand() {
					return StrandedFeature.UNKNOWN;
				}
			});
			
			i++;
		}
		writer.endDocument();

	}
}