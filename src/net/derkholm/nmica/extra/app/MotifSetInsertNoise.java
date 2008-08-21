package net.derkholm.nmica.extra.app;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.metamotif.Dirichlet;
import net.derkholm.nmica.model.motif.NMWeightMatrix;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.seq.NMSimpleDistribution;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Add noise to columns of a motif or remove/replace columns", 
		generateStub=true)
@NMExtraApp(launchName="nmnoise", vm=VirtualMachine.SERVER)
public class MotifSetInsertNoise {

	private double maxPerturb = 0.5;
	private int removeColumnsFromEnds = 0;
	private String outFilename;
	private int randomCols = 0;

	@Option(help="The maximum perturbation added to a nucleotide weight", optional=true)
	public void setMaxRate(double d) {
		this.maxPerturb = d;
	}
	
	@Option(help="Remove specified number of columns from ends randomly", optional=true)
	public void setRemoveCols(int i) {
		this.removeColumnsFromEnds = i;
	}
	
	@Option(help="Sample a random column", optional=true)
	public void setRandomCols(int i) {
		this.randomCols = i;
	}
	
	@Option(help="Output filename", optional=true)
	public void setOut(String str) {
		this.outFilename = str;
	}
	
	private Motif[] motifs;
	
	private static Random generator = new Random();
	
	public void main(String[] args) throws FileNotFoundException, Exception {
		List<Motif> motifList = new ArrayList<Motif>();
		for (String filen : args) {
			Motif[] ms = MotifIOTools.loadMotifSetXML(
				new BufferedInputStream(new FileInputStream(filen)));
			
			for (Motif m : ms) motifList.add(m);
		}
		
		this.motifs = motifList.toArray(new Motif[motifList.size()]);
		
		List<Motif> outputMotifs = new ArrayList<Motif>();
		
		for (Motif m : this.motifs) {
			Motif nm = new Motif();
			nm.setName(m.getName());
			nm.setThreshold(m.getThreshold());
			nm.setVersion(m.getVersion());
			
			WeightMatrix wm;
			wm = addNoise(m.getWeightMatrix(), this.maxPerturb);
			
			if (removeColumnsFromEnds > 0) {
				for (int i = 0; i < removeColumnsFromEnds; i++) {
					wm = removeColumnFromEitherEnd(wm);
				}
			}
			
			if (randomCols > 0) {
				if (randomCols < wm.columns()) {
					System.err.printf(
							"randomCols > motif length (%d > %d)%n",
							randomCols,wm.columns());
					System.exit(1);
				}
				Dirichlet dir = new Dirichlet(DNATools.getDNA());
				
				Set<Integer> replacedColumns = new TreeSet<Integer>();
				
				for (int i = 0; i < randomCols; i++) {
					int repCol = -1;
					do {
						repCol = generator.nextInt(wm.columns());
					} while (replacedColumns.contains(repCol));
					
					replacedColumns.add(repCol);
					
					wm = replaceColumnWithSampleFromDirichlet(
							wm,repCol,dir);
				}
			}
			
			nm.setWeightMatrix(wm);
			outputMotifs.add(nm);
		}
		
		if (outFilename != null)
			MotifIOTools.writeMotifSetXML(new BufferedOutputStream(
											new FileOutputStream(outFilename)), 
											outputMotifs.toArray(
													new Motif[outputMotifs.size()]));
		else {
			MotifIOTools.writeMotifSetXML(System.out, 
					outputMotifs.toArray(
							new Motif[outputMotifs.size()]));
		}
	}
	
	private WeightMatrix removeColumnFromEitherEnd(WeightMatrix wm) {
		Distribution[] ds = new Distribution[wm.columns() - 1];
		if (Math.random() < 0.5) {
			for (int i = 0; i < wm.columns()-1; i++) {
				ds[i] = new SimpleDistribution(wm.getColumn(i));
			}
		} else {
			for (int i = 1; i < wm.columns(); i++) {
				ds[i] = new SimpleDistribution(wm.getColumn(i));
			}
		}
		try {
			return new NMWeightMatrix(ds,ds.length,0);
		} catch (IllegalAlphabetException e) {
			throw new BioError("Illegal alphabet exception caught",e);
		}
	}
	
	private WeightMatrix replaceColumnWithSampleFromDirichlet(
			WeightMatrix wm, int col, Dirichlet dir) {
		Distribution[] ds = new Distribution[wm.columns()];
		
		for (int i = 0; i < ds.length; i++)
			if (i == col) ds[i] = dir.sampleDistribution();
			else ds[i] = new SimpleDistribution(ds[i]);
		
		try {
			return new NMWeightMatrix(ds,ds.length,0);
		} catch (IllegalAlphabetException e) {
			throw new BioError("Illegal alphabet exception caught",e);
		}
	}
	
	private WeightMatrix addNoise(WeightMatrix wm, double maxRate) throws IllegalAlphabetException {
		Distribution[] ds = new Distribution[wm.columns()];
		for (int i = 0; i < wm.columns(); i++) {
			Distribution d = wm.getColumn(i);
			Distribution noisyD = addNoise(d, maxRate);
			
			ds[i] = noisyD;
		}
		
		WeightMatrix nwm = new NMWeightMatrix(ds, ds.length, 0);
		return nwm;
	}
	
	private Distribution addNoise(Distribution old, double maxRate) {
		try {
            if (Math.random() < 0.5) {
                return sampleDist(old, maxRate);
            } else {
                return sampleDistBack(old, maxRate);
            }
        } catch (Exception ex) {
            throw new BioError(ex);
        }
	}
	
	private Distribution sampleDist(Distribution dist, double step) 
       throws Exception {
       FiniteAlphabet alpha = (FiniteAlphabet) dist.getAlphabet();
       Distribution nudist = new NMSimpleDistribution(alpha);
       AtomicSymbol beneficiary = (AtomicSymbol) new UniformDistribution(alpha).sampleSymbol();
       double inc = dist.getWeight(beneficiary) * step;
       for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
           AtomicSymbol s = (AtomicSymbol) si.next();
           double p = dist.getWeight(s);
           double w;
           if (s == beneficiary) {
               w = (p + inc) / (1.0 + inc);
           } else {
               w = p / (1.0 + inc);
           }
           nudist.setWeight(s, w);
       }
       return nudist;
   }
   
   private Distribution sampleDistBack(Distribution dist, double step) 
       throws Exception
   {
       FiniteAlphabet alpha = (FiniteAlphabet) dist.getAlphabet();
       Distribution nudist = new NMSimpleDistribution(alpha);
       AtomicSymbol beneficiary = (AtomicSymbol) new UniformDistribution(alpha).sampleSymbol();
       
       double bFwd = dist.getWeight(beneficiary);
       double bBack = bFwd / (1.0 + step * (1 - bFwd));
       for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
           AtomicSymbol s = (AtomicSymbol) si.next();
           double p = dist.getWeight(s);
           double w;
           if (s == beneficiary) {
               w = bBack;
           } else {
               w = p * (1.0 - bBack) / (1.0 - bFwd);
           }
           nudist.setWeight(s, w);
       }
       return nudist;
   }
}