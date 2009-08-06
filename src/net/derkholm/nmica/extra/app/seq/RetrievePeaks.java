package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import biobits.utils.IOTools;

@App(overview = "Get sequences from Ensembl that using peaks in the FindPeaks format", generateStub = true)
@NMExtraApp(launchName = "nmpeakseq", vm = VirtualMachine.SERVER)
public class RetrievePeaks extends RetrieveEnsemblSequences {

	private File peaksFile;
	private File outFile;

	@Option(help="Peaks")
	public void setPeaks(File f) {
		peaksFile = f;
	}
	
	@Option(help="Output file")
	public void setOut(File f) {
		this.outFile = f;
	}
	
	public static class PeakEntry implements Comparable<PeakEntry> {
		public final int id;
		public String seqName;
		public final int startCoord;
		public final int endCoord;
		public final int regEndCoord;
		public final double value;
		
		public PeakEntry(
				int id, 
				String seqName,
				int startCoord, 
				int endCoord, 
				int regEndCoord, 
				double value) {
			this.id = id;
			this.seqName = seqName;
			this.startCoord = startCoord;
			this.endCoord = endCoord;
			this.regEndCoord = regEndCoord;
			this.value = value;
		}

		public int compareTo(PeakEntry o) {
			return Double.compare(this.value, o.value);
		}
	}
	
	public void main(String[] argv) throws Exception {
		initializeEnsemblConnection();
		
		BufferedReader br = IOTools.inputBufferedReader(argv);
		
		Map<String,List<Location>> locations = new HashMap<String,List<Location>>();
		String line = null;
		List<PeakEntry> peaks = new ArrayList<PeakEntry>();
		while ((line = br.readLine()) != null) {
			StringTokenizer tok = new StringTokenizer("\t");
			
			int id = Integer.parseInt(tok.nextToken());
			String chromo = tok.nextToken();
			int startCoord = Integer.parseInt(tok.nextToken());
			int endCoord = Integer.parseInt(tok.nextToken());
			int regEndCoord = Integer.parseInt(tok.nextToken());
			double value = Double.parseDouble(tok.nextToken());
			
			peaks.add(new PeakEntry(id, chromo, startCoord, endCoord, regEndCoord, value));
		}
		
		Map<String, List<Location>> locationMap = new HashMap<String,List<Location>>();
		
		for (PeakEntry peak : peaks) {
			if (locationMap.get(peak.seqName) == null) {
				
			}
			locationMap.get(peak.seqName).add(new RangeLocation(peak.startCoord,peak.endCoord));
		}
		
		for (String chromo : locationMap.keySet()) {
			Location loc = LocationTools.union(locationMap.get(chromo));
			
			for (Iterator<?> bi = loc.blockIterator(); bi.hasNext();) {
				RangeLocation l = (RangeLocation)bi.next();
				
				FeatureHolder feats = 
					seqDB.filter(
						new FeatureFilter.ContainedByLocation(l));
				
				for (Iterator<?> fi = feats.features(); fi.hasNext();) {
					StrandedFeature peakFeat = (StrandedFeature) fi.next();
					Sequence chr = peakFeat.getSequence();
					
					SymbolList symList = chr.subList(l.getMin(), l.getMax());
					Sequence seq = new SimpleSequence(symList, null, 
							String.format("%s_%d-%d", 
									chromo, 
									l.getMin(),
									l.getMax()),
							Annotation.EMPTY_ANNOTATION);
					
					RichSequence.IOTools.writeFasta(System.out, seq, null);
				}
			}
		}
	}
}