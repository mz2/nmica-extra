package net.derkholm.nmica.extra.app.motif;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.motif.SquaredDifferenceMotifComparitor;
import net.derkholm.nmica.motif.align.MotifAlignment;
import net.derkholm.nmica.motif.align.MotifPair;
import net.derkholm.nmica.seq.WmTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@SuppressWarnings("unchecked")
@App(overview="Assess significance of motif comparison", generateStub=true)
@NMExtraApp(launchName="nmshuffle", vm=VirtualMachine.SERVER)
public class MotifSetComparisonShuffler {

	private int bootstraps = 100;
	private Shuffler shuffle = Shuffler.SHUFFLE;
	private boolean ignoreDiagonal = false;
	private double threshold = Double.POSITIVE_INFINITY;
	private File database;
	private File motifs;
	private boolean outputPairedMotifs;
	private File outputFile;
	
	@Option(help="Significance threshold for reporting matches", optional=true)
	public void setThreshold(double d) {
		this.threshold = d;
	}
	
	@Option(help="Ignore target motifs with matching names", optional=true)
	public void setIgnoreDiagonal(boolean b) {
		this.ignoreDiagonal = b;
	}
	
	@Option(help="Shuffling method", optional=true)
	public void setShuffle(Shuffler shuffle) {
		this.shuffle = shuffle;
	}
	
	@Option(help="Number of shufflings", optional=true)
	public void setBootstraps(int i) {
		this.bootstraps = i;
	}
	
	@Option(help="Database of motifs to compare motif set against (can also alternatively supply database and query motifs as unnamed arguments)", optional=true)
	public void setDatabase(File f) {
		this.database = f;
	}
	
	@Option(help="Source motif set (can also alternatively supply database and query motifs as unnamed arguments)", optional=true)
	public void setMotifs(File f) {
		this.motifs = f;
	}
	
	@Option(help="Output paired motifs", optional=true)
	public void setOutputPairedMotifs(boolean b) {
		this.outputPairedMotifs = b;
	}
	
	@Option(help="Output file (written to stdout if blank)", optional=true)
	public void setOut(File f) {
		this.outputFile = f;
	}
	
	/**
	 * @param args
	 */
	public void main(String[] args) 
		throws Exception
	{
		Motif[] dictionary, queries;
		if (args.length == 2) {
			dictionary = MotifIOTools.loadMotifSetXML(new FileReader(args[0]));
			queries = MotifIOTools.loadMotifSetXML(new FileReader(args[1]));
		} else {
			dictionary = MotifIOTools.loadMotifSetXML(new FileReader(this.database));
			queries = MotifIOTools.loadMotifSetXML(new FileReader(this.motifs));
		}
		
		PrintStream outputStream = null;

		if (outputFile == null) {
			outputStream = System.out;
		} else {
			outputStream = new PrintStream(new FileOutputStream(outputFile));
		}
		
		Distribution elsewhere = new UniformDistribution(DNATools.getDNA());
		
		List<MotifPair> motifPairs = new ArrayList<MotifPair>();
		
		for (Motif query : queries) {
			System.err.println("Searching " + query.getName());
			double bestScore = Double.POSITIVE_INFINITY;
			Motif bestMatch = null;
			for (Motif m : dictionary) {
				if (ignoreDiagonal && m.getName().equals(query.getName())) {
					continue;
				}
				double score = diCompareMotifs(query.getWeightMatrix(), elsewhere, m.getWeightMatrix(), elsewhere);
				if (score < bestScore) {
					bestScore = score;
					bestMatch = m;
				}
			}
			
			double[] bsScores = new double[bootstraps];
			for (int b = 0; b < bootstraps; ++b) {
				WeightMatrix sq = shuffle.shuffle(query.getWeightMatrix());
				double mbs = Double.POSITIVE_INFINITY;
				for (Motif m : dictionary) {
					if (ignoreDiagonal && m.getName().equals(query.getName())) {
						continue;
					}
					mbs = Math.min(mbs, diCompareMotifs(sq, elsewhere, m.getWeightMatrix(), elsewhere));
				}
				bsScores[b] = mbs;
			}
			Arrays.sort(bsScores);
			
			int c = 0;
			while (c < bootstraps && bsScores[c] < bestScore) {
				++c;
			}
			
			double pv = ((1.0 * c) / bootstraps);
			if (pv <= threshold) {
				if (!outputPairedMotifs) {
					outputStream.printf("%s\t%s\t%g\t%d\t%g%n", query.getName(), bestMatch.getName(), bestScore, c, pv);
				} else {
					motifPairs.add(new MotifPair(query,bestMatch,bestScore,pv,false));
					
					bestMatch.getAnnotation().setProperty("reciprocal_match", query.getName());
					bestMatch.getAnnotation().setProperty("reciprocal_match_p-value", pv);
					bestMatch.getAnnotation().setProperty("reciprocal_match_best_score", bestScore);

					query.getAnnotation().setProperty("reciprocal_match", bestMatch.getName());
					query.getAnnotation().setProperty("reciprocal_match_p-value", pv);
					query.getAnnotation().setProperty("reciprocal_match_best_score", bestScore);
				}
			}
		}
		
		if (motifPairs.size() > 0) {
			List<Motif> motifs = new ArrayList<Motif>();
			
			for (MotifPair mp : motifPairs) {
				MotifAlignment alignment 
					= new MotifAlignment(
						new Motif[]{mp.getM1(),mp.getM2()}, 
						SquaredDifferenceMotifComparitor.getMotifComparitor());
				
				alignment = new MotifAlignment(alignment.motifs(), SquaredDifferenceMotifComparitor.getMotifComparitor());
				alignment = new MotifAlignment(alignment.motifs(), SquaredDifferenceMotifComparitor.getMotifComparitor());
				
				motifs.add(alignment.motifs()[0]);
				motifs.add(alignment.motifs()[1]);
			}
			
			MotifIOTools.writeMotifSetXML(
					outputStream, 
					motifs.toArray(new Motif[motifs.size()]));
		}
	}
	
	private final static List<Symbol> symL;
	
	static {
		try {
			symL = (List<Symbol>) DNATools.createDNA("acgt").toList();
		} catch (IllegalSymbolException e) {
			throw new Error(e);
		}
	}
	
	private static WeightMatrix shuffleWM(WeightMatrix oldm)
		throws Exception
	{
		List l = new ArrayList();
        for (int c = 0; c < oldm.columns(); ++c) {
            l.add(oldm.getColumn(c));
        }
        Collections.shuffle(l);
        return new SimpleWeightMatrix((Distribution[]) l.toArray(new Distribution[0]));
	}
	
	private static WeightMatrix randomizeWM(WeightMatrix oldm)
		throws Exception
	{
	    Distribution[] dists = new Distribution[oldm.columns()];
	    for (int c = 0; c < oldm.columns(); ++c) {
	        Distribution oDist = oldm.getColumn(c);
	        Distribution nDist = DistributionFactory.DEFAULT.createDistribution(oDist.getAlphabet());
	        List newL =  new ArrayList(symL);
	        Collections.shuffle(newL);
	        for (int i = 0; i < symL.size(); ++i) {
	            nDist.setWeight((Symbol) newL.get(i), oDist.getWeight((Symbol) symL.get(i)));
	        }
	        dists[c] = nDist;
	    }
	    return new SimpleWeightMatrix(dists);
	}
	
    private static double div(Distribution d0, Distribution d1) 
		throws Exception
    {
        double cScore = 0.0;
        for (Iterator i = ((FiniteAlphabet) d0.getAlphabet()).iterator(); i.hasNext(); ) {
             Symbol s= (Symbol) i.next();
             double delta = d0.getWeight(s) - d1.getWeight(s);
             cScore += (delta * delta);
        }
        return Math.pow(cScore, 2.5 / 2.0);
    }
    

    private static double diCompareMotifs(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1)
    	throws Exception
    {
    	return Math.min(
    			compareMotifs(wm0, pad0, wm1, pad1),
    			compareMotifs(wm0, pad0, WmTools.reverseComplement(wm1), pad1)
    	);
    }
    
    private static double compareMotifs(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1)
		throws Exception
	{
		double bestScore = Double.POSITIVE_INFINITY;
		int minPos = -wm1.columns();
		int maxPos = wm0.columns() + wm1.columns();
		for (int offset = -wm1.columns(); offset <= wm0.columns(); ++offset) {
			double score = 0.0;
			for (int pos = minPos; pos <= maxPos; ++pos) {
				Distribution col0 = pad0, col1 = pad1;
				if (pos >= 0 && pos < wm0.columns()) {
					col0 = wm0.getColumn(pos);
				}
				int opos = pos - offset;
				if (opos >= 0 && opos < wm1.columns()) {
					col1 = wm1.getColumn(opos);
				}
				double cScore = div(col0, col1);
                                score += cScore;
			}
			bestScore = Math.min(score, bestScore);
		}
		return bestScore;
	}

    private static enum Shuffler {
    	SHUFFLE {
    		public WeightMatrix shuffle(WeightMatrix wm) 
    			throws Exception
    		{
    			return shuffleWM(wm);
    		}
    	},
    	RANDOMIZE {
    		public WeightMatrix shuffle(WeightMatrix wm) 
    			throws Exception
    		{
    			return randomizeWM(wm);
    		}
    	},
    	BOTH {
    		public WeightMatrix shuffle(WeightMatrix wm) 
    			throws Exception
    		{
    			return randomizeWM(shuffleWM(wm));
    		}
    	};
    	
    	public abstract WeightMatrix shuffle(WeightMatrix wm) throws Exception;
    }
    
}
