package net.derkholm.nmica.extra.app;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifComparisonMatrixBundle;
import net.derkholm.nmica.motif.MotifComparitorIFace;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.motif.SimpleMotifComparitor;
import net.derkholm.nmica.motif.align.MotifAlignment;
import net.derkholm.nmica.motif.align.MotifPairWithOffset;

import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@App(overview="A tool for aligning a set of motifs and output the alignment", 
	generateStub=true)
@NMExtraApp(launchName="nmalign", vm=VirtualMachine.SERVER)
public class MotifAligner {
	
	private File motifs;

	@Option(help="Input motifset")
	public void setMotifs(File motifs) {
		this.motifs = motifs;
	}
	
	public void main(String[] args) throws FileNotFoundException, Exception {
		Motif[] motifs = MotifIOTools.loadMotifSetXML(new FileReader(this.motifs));
		List<Motif> motifList = new ArrayList<Motif>(Arrays.asList(motifs));
		
		SimpleMotifComparitor mc = SimpleMotifComparitor.getMotifComparitor();
		
		MotifAlignment alignment = new MotifAlignment((FiniteAlphabet)
					motifs[0].getWeightMatrix().getColumn(0).getAlphabet());
		
		//initialise by aligning the closest pair of motifs
		{ 
			MotifPairWithOffset mp = highestScoringPair((MotifComparitorIFace)mc, 
							motifList.toArray(new Motif[motifList.size()]));
			mp.getM1();
			mp.getM2();
			
			alignment.addMotif(mp.getM1(), 0);
			alignment.addMotif(mp.getM2(), mp.getOffset());
			
			System.out.println(	mp.getM1().getName() + " + " + 
								mp.getM2().getName() + 
								" flipped:" + mp.isFlipped());
			
			motifList.remove(mp.getM1());
			motifList.remove(mp.getM2());
		}
		
		//then iterate through the remaining ones
		while (motifList.size() > 0) {
			System.out.print("motifList:");
			for (Motif m : motifList)
				System.out.print(m.getName() + " ");
			System.out.println();
			
			System.out.print("aligned  :");
			for (Motif m : alignment.motifs())
				System.out.print(m.getName() + " ");
			System.out.println();
			
			Motif[] remainingMotifs 
				= motifList.toArray(new Motif[motifList.size()]);
			
			MotifPairWithOffset mp 
				= highestScoringPair(mc, alignment.motifs(), remainingMotifs);
			
			System.out.println(	mp.getM1().getName() + " - " + 
								mp.getM2().getName() + " flipped:" + mp.isFlipped());
			
			alignment.addMotif(
					mp.getM2(),
					alignment.offset(mp.getM1()) + mp.getOffset());
			
			motifList.remove(mp.getM2());
			System.out.println(motifList.size());
		}
	}
	
	public MotifPairWithOffset highestScoringPair(MotifComparitorIFace mc, Motif[] motifs) 
		throws IllegalAlphabetException, IllegalSymbolException {
		
		MotifComparisonMatrixBundle mb = mc.
			getComparisonMatrixWithOffsets(motifs);
	
		double maxScore = Double.POSITIVE_INFINITY;
		int offset = 0;
		Motif m0, m1;
		m0 = null;
		m1 = null;
		boolean flipped = false;
		
		for (int i = 0; i < motifs.length; i++) {
			for (int j = (i+1); j < motifs.length; j++) {
				double ijScore = mb.getSenseScoreMatrix().get(i, j);
				if (ijScore < maxScore) {
					maxScore = ijScore;
					offset = mb.getSenseOffsetMatrix().get(i, j);
					m0 = motifs[i];
					m1 = motifs[j];
				}
				
				double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
				if (fijScore < maxScore) {
					maxScore = fijScore;
					offset = mb.getAntisenseOffsetMatrix().get(i, j);
					m1 = motifs[j];
					m0 = motifs[i];
					flipped = true;
				}
			}
		}
		System.out.println("OFFSET:::"+offset);
		
		if (m0 != null && m1 != null) {
			MotifPairWithOffset mpoffset = new MotifPairWithOffset(
					m0, 
					m1, 
					maxScore, 
					flipped, offset);					
			return mpoffset;
		} else {
			throw new IllegalArgumentException("Could not align motifs");
		}
	}
	
	public MotifPairWithOffset highestScoringPair(
			MotifComparitorIFace mc, Motif[] motifs0, Motif[] motifs1) 
		throws Exception {
		MotifComparisonMatrixBundle mb = mc.
			getComparisonMatrixWithOffsets(
				motifs0, motifs1);
		
		double maxScore = Double.POSITIVE_INFINITY;
		int offset = 0;
		Motif m0, m1;
		m0 = null;
		m1 = null;
		boolean flipped = false;
		
		for (int i = 0; i < motifs0.length; i++) {
			for (int j = 0; j < motifs1.length; j++) {
				double ijScore = mb.getSenseScoreMatrix().get(i, j);
				if (ijScore < maxScore) {
					maxScore = ijScore;
					offset = mb.getSenseOffsetMatrix().get(i, j);
					m0 = motifs0[i];
					m1 = motifs1[j];
				}
				
				double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
				if (fijScore < maxScore) {
					maxScore = fijScore;
					offset = mb.getAntisenseOffsetMatrix().get(i, j);
					m1 = motifs0[j];
					m0 = motifs1[i];
					flipped = true;
				}
			}
		}
		
		if (m0 != null && m1 != null) {
			MotifPairWithOffset mpoffset = new MotifPairWithOffset(m0, m1, maxScore, flipped, offset);					
			return mpoffset;
		} else {
			throw new IllegalArgumentException("Could not align motifs");
		}

	}
}