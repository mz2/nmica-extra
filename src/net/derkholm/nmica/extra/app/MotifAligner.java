package net.derkholm.nmica.extra.app;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.extra.motif.comparison.BLiCMotifComparitor;
import net.derkholm.nmica.extra.motif.comparison.KullbackLeiblerDifferenceMotifComparitor;
import net.derkholm.nmica.model.metamotif.MetaMotif;
import net.derkholm.nmica.model.metamotif.MetaMotifIOTools;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifComparisonMatrixBundle;
import net.derkholm.nmica.motif.MotifComparitorIFace;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.motif.SquaredDifferenceMotifComparitor;
import net.derkholm.nmica.motif.align.MotifAlignment;
import net.derkholm.nmica.motif.align.MotifPairWithOffset;
import net.derkholm.nmica.motif.align.MotifAlignment.MotifAlignmentElement;
import net.derkholm.nmica.seq.WmTools;

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
	private String outputType = "align_cons";
	private String outFile;
	private MotifComparitorIFace mc = SquaredDifferenceMotifComparitor.getMotifComparitor();
	
	@Option(help="Input motifset")
	public void setMotifs(File motifs) {
		this.motifs = motifs;
	}
	
	@Option(help="Output file",optional=true)
	public void setOut(String outFile) {
		this.outFile = outFile;
	}
	
	@Option(help="Distance metric:sqdiff2.5(default)|sqdiff|kd|blic",optional=true)
	public void setDist(String dist) {
		if (!dist.equals("sqdiff2.5") &!
				dist.equals("sqdiff") &!
				dist.equals("kd") &!
				dist.equals("blic"))
			return;
		
		if (dist.equals("sqdiff2.5"))
			mc = SquaredDifferenceMotifComparitor.getMotifComparitor();
		else if (dist.equals("sqdiff"))
			mc = SquaredDifferenceMotifComparitor.getCarthesianMotifComparitor();
		else if (dist.equals("kd"))
			mc = KullbackLeiblerDifferenceMotifComparitor.getMotifComparitor();
		else if (dist.equals("blic"))
			mc = BLiCMotifComparitor.getMotifComparitor();
		else {
			System.err.println("Incorrect distance metric type : " + dist);
			System.exit(1);
		}
	}
	
	@Option(help="Output type (default=align_cons): avg|metamotif|cons|align_cons",optional=true)
	public void setOutputType(String type) {
		if (!type.equals("avg") &! 
				type.equals("all") &!
				type.equals("metamotif") &! 
				type.equals("cons") &!
				type.equals("align_cons")) {
			System.err.println("Allowed output types: avg|metamotif|cons|align_cons");
			System.exit(1);
		}
		this.outputType = type;
	}
	
	public void main(String[] args) throws FileNotFoundException, Exception {
		Motif[] motifs = MotifIOTools.loadMotifSetXML(new FileReader(this.motifs));
		List<Motif> motifList = new ArrayList<Motif>(Arrays.asList(motifs));

		MotifAlignment alignment = new MotifAlignment((FiniteAlphabet)
					motifs[0].getWeightMatrix().getColumn(0).getAlphabet());
		
		/*initialise by aligning the closest pair of motifs*/
		{ 
			MotifPairWithOffset mp = highestScoringPair((MotifComparitorIFace)mc, 
							motifList.toArray(new Motif[motifList.size()]));
			mp.getM1();
			mp.getM2();
			
			alignment.addMotif(mp.getM1(), 0, false);
			alignment.addMotif(mp.getM2(), mp.getOffset(), mp.isFlipped());
			
			/*if (mp.isFlipped())
				alignment.addMotif(mp.getM2(), -mp.getOffset(), true);
			else
				alignment.addMotif(mp.getM2(), mp.getOffset(), false);
			*/
			
			/* if it's flipped, let's flip it back */
			if (mp.isFlipped()) {
				//System.err.println("|FLIPPED " + mp.getM2().getName() + "|");
				mp.getM2()
					.setWeightMatrix(
						WmTools.reverseComplement(
							mp.getM2().getWeightMatrix()));
			}
			
			/*
			System.err.println(	" | " + mp.getM1().getName() + " - " + 
								mp.getM2().getName() + 
								" flipped:" + mp.isFlipped() + 
								" offset: " + mp.getOffset() + 
								" score:" + mp.getScore());
			*/
			
			motifList.remove(mp.getM1());
			motifList.remove(mp.getM2());
		}
		
		/* then iterate through the remaining ones */
		while (motifList.size() > 0) {
			/*
			System.err.print("motifList:");
			for (Motif m : motifList)
				System.err.print(m.getName() + " ");
			System.err.println();
			
			System.err.print("aligned  :");
			for (Motif m : alignment.motifs())
				System.err.print(m.getName() + " ");
			System.err.println();
			*/
			
			Motif[] remainingMotifs 
				= motifList.toArray(new Motif[motifList.size()]);
			
			MotifPairWithOffset mp 
				= highestScoringPair(mc, alignment.motifs(), remainingMotifs);
			
			/* if it's flipped, let's flip it back */
			if (mp.isFlipped()) {
				//System.err.println("**FLIPPED " + mp.getM2().getName() + "**");
				mp.getM2()
					.setWeightMatrix(
						WmTools.reverseComplement(
							mp.getM2().getWeightMatrix()));
			}
			
			/*
			System.err.println(	" > " + mp.getM1().getName() + " - " + 
								mp.getM2().getName() + 
								" flipped:" + mp.isFlipped() + 
								" offset: " + mp.getOffset() + 
								" score:" + mp.getScore());
			*/
			int offset;
			if (mp.isFlipped()) {
				offset = alignment.offset(mp.getM1()) + mp.getOffset();
			} else {
				offset = alignment.offset(mp.getM1()) + mp.getOffset();
			}
			/*
			System.err.println("m1 offset:" + alignment.offset(mp.getM1()));
			System.err.println("mp offset:" + mp.getOffset());
			System.err.println("OFFSET:" + offset);
			*/
			
			alignment.addMotif(
						mp.getM2(),
						offset,
						mp.isFlipped());
			
			motifList.remove(mp.getM2());
		}

		/*
		for (MotifAlignmentElement elem : alignment)
			System.err.println(elem.getMotif().getName() + " " + elem.getOffset() + " " + elem.isFlipped());
		*/
		
		alignment = alignment.alignmentWithZeroOffset();
		
		if (outputType.equals("avg"))
			MotifIOTools.writeMotifSetXML(System.out, new Motif[] {alignment.averageMotif()});
		if (outputType.equals("all"))
			MotifIOTools.writeMotifSetXML(System.out, alignment.motifs());
		else if (outputType.equals("metamotif"))
			MetaMotifIOTools.writeMetaMotifSetToMotifSetWithAnnotations(System.out, new MetaMotif[] {alignment.metamotif()});
		else if (outputType.equals("cons"))
			System.out.println(alignment.consensusString());
		else if (outputType.equals("align_cons"))
			System.out.println(alignment.alignmentConsensusString());

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
					flipped = false;
				}
				
				double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
				if (fijScore < maxScore) {
					maxScore = fijScore;
					offset = mb.getAntisenseOffsetMatrix().get(i, j);
					m0 = motifs[i];
					m1 = motifs[j];
					flipped = true;
				}
			}
		}
		//System.out.println("OFFSET:::"+offset);
		
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
		
		//System.err.println("---");
		for (int i = 0; i < motifs0.length; i++) {
			for (int j = 0; j < motifs1.length; j++) {
				double ijScore = mb.getSenseScoreMatrix().get(i, j);
				if (ijScore < maxScore) {
					maxScore = ijScore;
					offset = mb.getSenseOffsetMatrix().get(i, j);
					//System.err.println("score::"+ijScore + " offset:" + offset);
					m0 = motifs0[i];
					m1 = motifs1[j];
					flipped = false;
				}
				
				double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
				if (fijScore < maxScore) {
					maxScore = fijScore;
					offset = mb.getAntisenseOffsetMatrix().get(i, j);
					//System.err.println("f score::"+ijScore + " offset:" + offset);
					m0 = motifs0[i];
					m1 = motifs1[j];
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