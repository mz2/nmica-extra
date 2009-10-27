package net.derkholm.nmica.extra.app.motif;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.extra.motif.NotEnoughColumnsLeftException;
import net.derkholm.nmica.model.motif.NMWeightMatrix;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Removes columns from motifs that are below a certain information content", generateStub=true)
@NMExtraApp(launchName = "nmentropycutoff")
public class MotifSetInformationContentCutoff {

	private File outFile;
	private double cutoff;

	@Option(help="Output file (XMS). Outputs to stdout if not specified.", optional=true)
	public void setOut(File f) {
		this.outFile = f;
	}
	
	@Option(help="Information content cutoff (default = 0.1)")
	public void setCutoff(double cutoff) {
		this.cutoff = cutoff;
	}
	
	public void main(String[] args) throws FileNotFoundException, Exception {

		List<Motif> outputMotifs = new ArrayList<Motif>();
		for (String arg : args) {
			Motif[] motifs = MotifIOTools.loadMotifSetXML(new FileReader(new File(arg)));
			
			for (Motif m : motifs) {
				try {
					outputMotifs.add(filterColumnsByInformationContent(m));					
				} catch (NotEnoughColumnsLeftException exc) {
					System.err.println(exc.getMessage());
				}
			}
		}
		
		OutputStream outputStream = System.out;
		if (outFile != null) {
			outputStream = new FileOutputStream(this.outFile);
		}
		MotifIOTools.writeMotifSetXML(
			outputStream, 
				outputMotifs.toArray(
					new Motif[outputMotifs.size()]));
		outputStream.flush();
		outputStream.close();
	}
	
	public Motif filterColumnsByInformationContent(Motif m) 
		throws IllegalAlphabetException, NotEnoughColumnsLeftException {
		
		List<Distribution> dists = new ArrayList<Distribution>();
		
		Motif newM = new Motif(m);
		
		for (int i = 0; i < newM.getWeightMatrix().columns(); i++) {
			Distribution d = newM.getWeightMatrix().getColumn(i);
			dists.add(d);
		}
		
		Distribution elsewhere = new UniformDistribution((FiniteAlphabet)m.getWeightMatrix().getAlphabet());
		double entropyElsewhere = DistributionTools.totalEntropy(elsewhere);
		
		//remove columns from left
		int culledFromLeft = 0;
		while (dists.size() > 0) {
			Distribution d = dists.get(0);
			if ((entropyElsewhere - DistributionTools.totalEntropy(d)) < cutoff) {
				dists.remove(0);
			} else {
				break;
			}
			culledFromLeft += 1;
		}
		//remove columns from right
		int culledFromRight = 0;
		while (dists.size() > 0) {
			Distribution d = dists.get(dists.size()-1);
			if ((entropyElsewhere - DistributionTools.totalEntropy(d)) < cutoff) {
				dists.remove(dists.size()-1);
			} else {
				break;
			}
			culledFromRight += 1;
		}
		
		if (dists.size() < 2) {
			throw new NotEnoughColumnsLeftException(
				String.format(
						"motif %s has %d columns after culling with information content cutoff %.3f.",
						m.getName(),
						m.getWeightMatrix().columns(),
						this.cutoff));
		}
		
		newM.setWeightMatrix(
				new NMWeightMatrix(
					(Distribution[]) dists.toArray(
						new Distribution[dists.size()]), dists.size(), 0));
		return newM;
	}
}
