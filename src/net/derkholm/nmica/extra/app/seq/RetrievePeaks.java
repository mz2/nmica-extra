package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.symbol.Location;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import biobits.utils.IOTools;

@App(overview = "Get sequences from Ensembl that using peaks in the FindPeaks format", generateStub = true)
@NMExtraApp(launchName = "nmpeakseq", vm = VirtualMachine.SERVER)
public class RetrievePeaks extends RetrieveEnsemblSequences {

	private File peaksFile;

	@Option(help="Peaks")
	public void setPeaks(File f) {
		peaksFile = f;
	}
	
	private static class PeakEntry implements Comparable<PeakEntry> {
		public final int id;
		public final int startCoord;
		public final int endCoord;
		public final int regEndCoord;
		public final double value;
		
		public PeakEntry(
				int id, 
				int startCoord, 
				int endCoord, 
				int regEndCoord, 
				double value) {
			this.id = id;
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
			
			peaks.add(new PeakEntry(id, startCoord, endCoord, regEndCoord, value));
		}
		
	}
}