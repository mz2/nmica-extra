package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.model.motif.extra.ScoredHit;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.seq.WmTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

//m.getName(), seq.getName(), max, (1.0 * gte) / bootstraps, Math.log10((1.0 * gte) / bootstraps));
@App(overview="Empirically find E-values for maxPerSeq motif scores.\nOutput:motif seq maxscore e-value log10(e-val)", generateStub=true)
@NMExtraApp(launchName = "nmempeval")
public class MotifSetEmpiricalEValueCalculator {

	private static final double LOG_2 = Math.log(2.0);
	
	private Motif[] motifs;
	private File seqs;
	private int bootstraps = 10000;
	private double pthresh = 1.0;

	private boolean collectHits;

	private ArrayList<ScoredHit> collectedHits;

	private boolean positiveHits;
	
	@Option(help="Number of bootstraps", optional=true)
	public void setBootstraps(int bootstraps) {
		this.bootstraps = bootstraps;
	}
	
	@Option(help="Probability threshold (default = 1.0)", optional=true)
	public void setThreshold(double d) {
		this.pthresh = d;
	}

	@Option(help="Motif set file")
	public void setMotifs(Reader min) 
		throws Exception
	{
		this.motifs = MotifIOTools.loadMotifSetXML(min);
	}

	@Option(help="Sequences to score", optional=true)
	public void setSeqs(File seqs) {
		this.seqs = seqs;
	}
	
	public void setCollectHits(boolean b) {
		this.collectHits = b;
		if (collectHits) {
			if (this.collectedHits == null) {
				this.collectedHits = new ArrayList<ScoredHit>();
			}
		}
	}
	
	public List<ScoredHit> collectedHits() {
		return this.collectedHits;
	}
	
	public void clearCollectedHits() {
		this.collectedHits = new ArrayList<ScoredHit>();
	}
	
	/* Only used when collectHits = true. Assigns a label to the hits 
	 * (whether they're from the positive or the negative set) */
	public void setPositiveHits(boolean b) {
		this.positiveHits = b;
	}

	public void calculate() throws Exception {
		AlphabetIndex index = AlphabetManager.getAlphabetIndex(DNATools.getDNA());
		
		for (Motif m : motifs) {
			Scanner fScanner = makeScanner(index, m.getWeightMatrix());
			Scanner rScanner = makeScanner(index, WmTools.reverseComplement(m.getWeightMatrix()));
			for (SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(seqs))); si.hasNext(); ) {
				Sequence seq = si.nextSequence();
				byte[] sin = indexSeq(index, seq);
				
				double max = Math.max(
						maxScore(fScanner, sin),
						maxScore(rScanner, sin)
				);
				
				byte[] rSin = new byte[sin.length];
				System.arraycopy(sin, 0, rSin, 0, sin.length);
				
				int gte = 0;
				int gteThresh = (int) Math.ceil(pthresh * bootstraps);
				for (int rep = 0; rep < bootstraps && gte < gteThresh; ++rep) {
					shuffle(rSin);
					if (isMaxScoreGreater(fScanner, rSin, max) || isMaxScoreGreater(rScanner, rSin, max)) {
						++gte;
					}
				}
				
				if (!collectHits) {
					System.out.printf("%s\t%s\t%g\t%g\t%g%n", 
							m.getName(), 
							seq.getName(), 
							max, 
							(1.0 * gte) / bootstraps, 
							Math.log10((1.0 * gte) / bootstraps));					
				} else {
					collectedHits.add(
						new ScoredHit(
							m.getName(),
							seq.getName(),
							positiveHits,
							max,
							(1.0 * gte) / bootstraps));
				}
			}
		}
	}
	/**
	 * @param args
	 */
	public void main(String[] args) 
		throws Exception
	{
		this.calculate();
	}

	private double maxScore(Scanner s, byte[] sin)
		throws Exception
	{
		double max = Double.NEGATIVE_INFINITY;
		int maxPos = sin.length - s.length();
		for (int p = 0; p < maxPos; ++p) {
			max = Math.max(max, s.score(sin, p, false));
		}
		double nm = max - s.getMaxScore();
		return nm;
	}
	
	private boolean isMaxScoreGreater(Scanner s, byte[] sin, double target)
		throws Exception
	{
		double max = Double.NEGATIVE_INFINITY;
		int maxPos = sin.length - s.length();
		for (int p = 0; p < maxPos; ++p) {
			if(s.score(sin, p, false) - s.getMaxScore() >= target) {
				return true;
			}
		}
		return false;
	}
	
	Random r = new Random();
	private void shuffle(byte[] ba) {
		int len = ba.length;
		for (int c = len; c > 1; --c) {
			int alt = r.nextInt(c);
			byte tmp = ba[c - 1];
			ba[c - 1] = ba[alt];
			ba[alt] = tmp;
		}
	}
	
	
    private static class Scanner {
        private Matrix2D bm;
        private AlphabetIndex index;
        private String name;
        private double maxScore;
        private int endPos;
        
        public int length() {
        	return bm.columns();
        }
        
        public int endPos() {
            return endPos;
        }
        
        Scanner(Matrix2D bm, AlphabetIndex index, String name) 
            throws Exception
        {
            this.bm = bm;
            this.index = index;
            this.name = name;
            this.maxScore = bmMaxScore(bm);
        }
        
        public double getMaxScore() {
            return maxScore;
        }
        
        public double score(byte[] sl, int pos, boolean gm)
            throws Exception
        {
            return scoreBM(sl, bm, pos, gm);
        }
        
        private  double scoreBM(byte[] sl, Matrix2D bm, int pos, boolean gm)
                throws Exception
        {
            double score = 0.0;
            int scol = 0;
            try {
                for (int col = 0; col < bm.columns(); ++col, ++scol) {
                    byte s = sl[pos + scol];
                    while (gm && scol > 0 && s < 0) {
                        s = sl[pos + ++scol];
                    }
                    if (s >= 0) {
                        score += bm.get(s, col);
                    } else {
                        return Double.NEGATIVE_INFINITY;
                    }
                }
            } catch (IndexOutOfBoundsException ex) {
                return Double.NEGATIVE_INFINITY;
            }
            endPos = pos + scol - 1;
            return score;
        }
        
        
        private double bmMaxScore(Matrix2D bm)
                throws Exception
        {
            double wmScore = 0.0;
            for (int c = 0; c < bm.columns(); ++c) {
                double colScore = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < bm.rows(); ++i) {
                    colScore = Math.max(colScore, bm.get(i, c));
                }
                wmScore += colScore;
            }
            return wmScore;
        }
    }
    
    private Scanner makeScanner(AlphabetIndex index, WeightMatrix wm)
    		throws Exception
    {
        FiniteAlphabet alpha = (FiniteAlphabet) wm.getAlphabet();
        Matrix2D bm = new SimpleMatrix2D(alpha.size(), wm.columns());
        for (int c = 0; c < wm.columns(); ++c) {
            Distribution wmCol = wm.getColumn(c);
            
            for (Iterator<?> si = alpha.iterator(); si.hasNext(); ) {
                Symbol s = (Symbol) si.next();
                double baseWeight = wmCol.getWeight(s);
                bm.set(index.indexForSymbol(s), c, Math.log(baseWeight) / LOG_2);
            }
        }
        return new Scanner(bm, index, null);
    }
    
    private byte[] indexSeq(AlphabetIndex index, SymbolList sl)
	    throws Exception
	{
	    Symbol gap = sl.getAlphabet().getGapSymbol();
	    
	    byte[] bsl = new byte[sl.length() + 1];
	    for (int i = 1; i <= sl.length(); ++i) {
	        Symbol s = sl.symbolAt(i);
	        if (s == gap) {
	            bsl[i] = -2;
	        } else if (s instanceof AtomicSymbol) {
	            bsl[i] = (byte) index.indexForSymbol(s);
	        } else {
	            bsl[i] = -1;
	        }
	    }
	    return bsl;
	}


}
