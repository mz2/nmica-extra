package net.derkholm.nmica.extra.app;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.derkholm.nmica.build.NMApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifAlignment;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.motif.SimpleMotifComparitor;
import net.derkholm.nmica.motif.SimpleMotifComparitor.MotifComparisonMatrixBundle;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@App(overview="A tool for aligning a set of motifs and output the alignment", 
	generateStub=true)
@NMApp(launchName="nmalign", vm=VirtualMachine.SERVER)
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
		
		while (motifList.size() > 0) {
			Motif[] remainingMotifs = motifList.toArray(new Motif[motifList.size()]);
			MotifComparisonMatrixBundle mb = mc.
												getComparisonMatrixWithOffsets(
													remainingMotifs);
			Motif m0 = null, m1 = null;
			double maxScore = Double.NEGATIVE_INFINITY;
			boolean flipped = false;
			int offset = 0;
			
			for (int i = 0; i < remainingMotifs.length; i++) {
				for (int j = (i+1); j < remainingMotifs.length; j++) {
					double ijScore = mb.getSenseScoreMatrix().get(i, j);
					if (ijScore > maxScore) {
						maxScore = ijScore;
						offset = mb.getSenseOffsetMatrix().get(i, j);
						m0 = remainingMotifs[i];
						m1 = remainingMotifs[j];
					}
					
					double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
					if (fijScore > maxScore) {
						maxScore = fijScore;
						offset = mb.getAntisenseOffsetMatrix().get(i, j);
						m0 = remainingMotifs[i];
						m1 = remainingMotifs[j];
						flipped = true;
					}
				}
			}
			
			if (m0 != null && m1 != null) {
				alignment.addMotif(m0, 0);
				alignment.addMotif(m1, offset);
			} else {
				throw new BioException("Could not align motifs");
			}
			
		}
	}
}