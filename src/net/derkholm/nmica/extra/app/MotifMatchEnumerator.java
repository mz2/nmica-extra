package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.maths.NativeMath;
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
	private WeightMatrix wm;
	private final static Symbol[] ALPHA = {DNATools.a(), DNATools.c(), DNATools.g(), DNATools.t()};
	
	@Option(help="The bit score threshold (default = -6.0)", optional=true)
	public void setScoreThreshold(double d) {
		this.scoreThreshold = d;
	}
	
	@Option(help="Input motif set (in the XMS format)", optional=true)
	public void setMotifs(File f) throws FileNotFoundException, Exception {
		Motif[] motifSet = MotifIOTools.loadMotifSetXML(new BufferedReader(new FileReader(f)));
		
		if (motifSet.length != 1) {
			System.err.printf("Motif set with exactly 1 member is expected as input.%n");
			System.exit(1);
		}
		
		wm = motifSet[0].getWeightMatrix();
	}
	
	/**
	 * @param args
	 */
	public void main(String[] args) throws Exception {
		rec(wm, new Symbol[wm.columns()], 0, 0);
	}

	private void rec(WeightMatrix wm, Symbol[] word, int index, double score)
		throws Exception
	{
		if (index == word.length) {
			StringBuilder sb = new StringBuilder();
			for (Symbol s : word) {
				sb.append(s.getName().charAt(0));
			}
			sb.append('\t');
			sb.append(score);
			System.out.println(sb.toString());
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
}
