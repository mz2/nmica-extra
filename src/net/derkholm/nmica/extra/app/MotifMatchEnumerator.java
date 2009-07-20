package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.analysis.ScoredString;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Compare real and theoretical motif-score distributions", generateStub=true)
@NMExtraApp(launchName="nmenummatches", vm=VirtualMachine.SERVER)
public class MotifMatchEnumerator {

	private double scoreThreshold = -6.0;
	private Motif motif;
	private boolean storeHits;
	private List<ScoredString> storedHits;
	private final static Symbol[] ALPHA = {DNATools.a(), DNATools.c(), DNATools.g(), DNATools.t()};
	
	@Option(help="The bit score threshold (default = -6.0)", optional=true)
	public void setScoreThreshold(double d) {
		this.scoreThreshold = d;
	}
	
	@Option(help="Input motif set in the XMS format. Needs to contain exactly one motif.", optional=true)
	public void setMotif(File f) throws FileNotFoundException, Exception {
		Motif[] motifSet = MotifIOTools.loadMotifSetXML(new BufferedReader(new FileReader(f)));
		
		if (motifSet.length != 1) {
			System.err.printf("Motif set with exactly 1 member is expected as input.%n");
			System.exit(1);
		}
		
		if (motifSet.length > 0) {
			System.err.printf("One motif expected in the input set. Got %d%n", motifSet.length);
			System.exit(1);
		}
		motif = motifSet[0];
	}
	
	@Option(help="Store results", optional=true)
	public void setStoreHits(boolean b) {
		this.storeHits = b;
		this.storedHits = new ArrayList<ScoredString>();
	}
	
	public List<ScoredString> storedHits() {
		return this.storedHits;
	}
	
	/**
	 * @param args
	 */
	public void main(String[] args) throws Exception {
		rec(motif.getWeightMatrix(), new Symbol[motif.getWeightMatrix().columns()], 0, 0);
	}

	private void rec(WeightMatrix wm, Symbol[] word, int index, double score)
		throws Exception
	{
		if (index == word.length) {
			StringBuilder sb = new StringBuilder();
			
			for (Symbol s : word) {sb.append(s.getName().charAt(0));}
			if (storeHits) {
				storedHits.add(new ScoredString(sb.toString(), score));
			} else {
				sb.append('\t');
				sb.append(score);
				System.err.println(sb.toString());
			}

		} else {
			Distribution d = wm.getColumn(index);
			double max = 0;
			for (Symbol s : ALPHA) {
				max = Math.max(max, d.getWeight(s));
			}
			for (Symbol s : ALPHA) {
				double here = score + Math.log(d.getWeight(s))/Math.log(2.0);
				if (here > scoreThreshold) {
					word[index] = s;
					rec(wm, word, index + 1, here);
				}
			}
		}
	}

	public void enumerateMatches(Motif m, double minThreshold) throws Exception {
		this.motif = m;
		this.setScoreThreshold(minThreshold);
		this.setStoreHits(true);
		this.main(null);
	}
}
