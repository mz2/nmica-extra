//
//  MotifComparitor.java
//  MotifExplorer
//
//  Created by Thomas Down on 08/03/2005.
//  Copyright 2005 __MyCompanyName__. All rights reserved.
//

package net.derkholm.nmica.motif.align;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import net.derkholm.nmica.matrix.IntMatrix2D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.matrix.SimpleIntMatrix2D;
import net.derkholm.nmica.matrix.SimpleMatrix2D;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.seq.WmTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;


public class SimpleMotifComparitor implements MotifComparitorIFace {
	private static SimpleMotifComparitor instance;
	
	public static SimpleMotifComparitor getMotifComparitor() {
		if (instance == null)
			instance = new SimpleMotifComparitor();
		
		return instance;
	}
	
    private Matrix2D[] getComparisonMatrix(Motif[] set0, Motif[] set1)
    	throws Exception {
        Matrix2D dMatrix = new SimpleMatrix2D(set0.length, set1.length);
        Matrix2D fMatrix = new SimpleMatrix2D(set0.length, set1.length);
        int numSteps = set0.length * set1.length;
        
        int step = 0;
        Distribution elsewhere = new UniformDistribution((FiniteAlphabet) set0[0].getWeightMatrix().getAlphabet());
        for (int i = 0; i < set0.length; ++i) {
            for (int j = 0; j < set1.length; ++j) {
                dMatrix.set(
                    i, j,
                    compareMotifs(set0[i].getWeightMatrix(), elsewhere, set1[j].getWeightMatrix(), elsewhere)
                );
                fMatrix.set(
                    i, j,
                    compareMotifs(set0[i].getWeightMatrix(), elsewhere, WmTools.reverseComplement(set1[j].getWeightMatrix()), elsewhere)
                );
                ++step;
            }
        }
        
        return new Matrix2D[] {dMatrix, fMatrix};
    }
    
    private Matrix2D[] getComparisonMatrix(Motif[] set) throws Exception {
        Matrix2D dMatrix = new SimpleMatrix2D(set.length, set.length);
        Matrix2D fMatrix = new SimpleMatrix2D(set.length, set.length);
        Distribution elsewhere = new UniformDistribution((FiniteAlphabet) set[0].getWeightMatrix().getAlphabet());
        
        for (int i = 0; i < set.length; ++i) {
            for (int j = (i+1); j < set.length; ++j) {
                dMatrix.set(
                    i, j,
                    compareMotifs(set[i].getWeightMatrix(), elsewhere, set[j].getWeightMatrix(), elsewhere)
                );
                fMatrix.set(
                    i, j,
                    compareMotifs(set[i].getWeightMatrix(), elsewhere, WmTools.reverseComplement(set[j].getWeightMatrix()), elsewhere)
                );
            }
        }
        
        return new Matrix2D[] {dMatrix, fMatrix};
    }
    
    public class MotifComparisonMatrixBundle {
    	public MotifComparisonMatrixBundle(Matrix2D matrix, Matrix2D matrix2,
				IntMatrix2D offsetMatrix, IntMatrix2D offsetMatrix2) {
			this.senseScoreMatrix = matrix;
			this.antisenseScoreMatrix = matrix2;
			this.senseOffsetMatrix = offsetMatrix;
			this.antisenseOffsetMatrix = offsetMatrix2;
		}
    	
		protected Matrix2D senseScoreMatrix;
    	protected Matrix2D antisenseScoreMatrix;
    	protected IntMatrix2D senseOffsetMatrix;
    	protected IntMatrix2D antisenseOffsetMatrix;
    	
		public Matrix2D getSenseScoreMatrix() {
			return senseScoreMatrix;
		}
		public void setSenseScoreMatrix(Matrix2D senseScoreMatrix) {
			this.senseScoreMatrix = senseScoreMatrix;
		}
		public Matrix2D getAntisenseScoreMatrix() {
			return antisenseScoreMatrix;
		}
		public void setAntisenseScoreMatrix(Matrix2D antisenseScoreMatrix) {
			this.antisenseScoreMatrix = antisenseScoreMatrix;
		}
		public IntMatrix2D getSenseOffsetMatrix() {
			return senseOffsetMatrix;
		}
		public void setSenseOffsetMatrix(IntMatrix2D senseOffsetMatrix) {
			this.senseOffsetMatrix = senseOffsetMatrix;
		}
		public IntMatrix2D getAntisenseOffsetMatrix() {
			return antisenseOffsetMatrix;
		}
		public void setAntisenseOffsetMatrix(IntMatrix2D antisenseOffsetMatrix) {
			this.antisenseOffsetMatrix = antisenseOffsetMatrix;
		}
    }
    
    public MotifComparisonMatrixBundle getComparisonMatrixWithOffsets(Motif[] set) throws Exception {
        Matrix2D dMatrix = new SimpleMatrix2D(set.length, set.length);
        Matrix2D fMatrix = new SimpleMatrix2D(set.length, set.length);
        IntMatrix2D dOffsetMatrix = new SimpleIntMatrix2D(set.length, set.length);
        IntMatrix2D fOffsetMatrix = new SimpleIntMatrix2D(set.length, set.length);
        
        Distribution elsewhere = new UniformDistribution((FiniteAlphabet) set[0].getWeightMatrix().getAlphabet());
        
        for (int i = 0; i < set.length; ++i) {
            for (int j = (i+1); j < set.length; ++j) {
            	ScoreOffsetPair dsp = compareMotifsWithOffset(set[i].getWeightMatrix(), elsewhere, set[j].getWeightMatrix(), elsewhere);
            	dMatrix.set(i, j, dsp.score);
            	dOffsetMatrix.set(i, j, dsp.offset);
            	
            	ScoreOffsetPair fsp = compareMotifsWithOffset(set[i].getWeightMatrix(), elsewhere, WmTools.reverseComplement(set[j].getWeightMatrix()), elsewhere);
            	fMatrix.set(i, j, fsp.score);
            	fOffsetMatrix.set(i, j, fsp.offset);
            }
        }
        
        return new MotifComparisonMatrixBundle(dMatrix, fMatrix,dOffsetMatrix,fOffsetMatrix);
    }

    private List<MotifPair> bestHits(Matrix2D dMatrix, Matrix2D fMatrix, Motif[] set0, Motif[] set1, boolean flip) {
         List<MotifPair> pairList = new ArrayList<MotifPair>();
         for (int i = 0; i < set0.length; ++i) {
	        int bestJ = -1;
	        double bestScore = Double.POSITIVE_INFINITY;
	        boolean bestIsFlipped = false;
	        for (int j = 0; j < set1.length; ++j) {
	            {
	                double score = dMatrix.get(i, j);
	                if (score < bestScore) {
	                    bestScore = score;
	                    bestJ = j;
	                    bestIsFlipped = false;
	                }
	            }
	            {
	                double score = fMatrix.get(i, j);
                    if (score < bestScore) {
	                    bestScore = score;
	                    bestJ = j;
	                    bestIsFlipped = true;
	                }
	            }
            }
            if (!flip) {
                pairList.add(new MotifPair(set0[i], set1[bestJ], bestScore, bestIsFlipped));
            } else {
                pairList.add(new MotifPair(set1[bestJ], set0[i], bestScore, bestIsFlipped));
            }
        }
        return pairList;
    }
    
    public Matrix2D bestHitsMatrix(Motif[] set0, Motif[] set1) 
    	throws Exception {
    	Matrix2D[] matrices = getComparisonMatrix(set0, set1);
    	Matrix2D dMatrix = matrices[0];
    	Matrix2D fMatrix = matrices[1];
    	
    	Matrix2D hitMatrix = new SimpleMatrix2D(set0.length, set1.length);
    	
    	for (int i = 0; i < set0.length; ++i) {
    		for (int j = 0; j < set1.length; j++) {
    			double d = dMatrix.get(i, j);
	    		double f = fMatrix.get(i, j);
	    		double best;
	    		boolean bestIsFlipped;
	    		if (d < f) {
	    			best = d; bestIsFlipped = false;
	    		} else {
	    			best = f; bestIsFlipped = true;
	    		}
	    		hitMatrix.set(i, j, best);
    		}
    	}
    	
    	return hitMatrix;
    }
    
    public Matrix2D bestHitsMatrix(Motif[] set) 
		throws Exception {
		Matrix2D[] matrices = getComparisonMatrix(set);
		Matrix2D dMatrix = matrices[0];
		Matrix2D fMatrix = matrices[1];
		
		Matrix2D hitMatrix = new SimpleMatrix2D(set.length, set.length);
		
		for (int i = 0; i < set.length; ++i) {
			for (int j = (i+1); j < set.length; j++) {
				double d = dMatrix.get(i, j);
	    		double f = fMatrix.get(i, j);
	    		double best;
	    		boolean bestIsFlipped;
	    		if (d < f) {
	    			best = d; bestIsFlipped = false;
	    		} else {
	    			best = f; bestIsFlipped = true;
	    		}
	    		hitMatrix.set(i, j, best);
			}
		}
		
		return hitMatrix;
	}
    
    public MotifPair[] allHits(Motif[] set0, Motif[] set1, double threshold)
	    throws Exception {
	    Matrix2D[] matrices = getComparisonMatrix(set0, set1);
	    Matrix2D dMatrix = matrices[0];
	    Matrix2D fMatrix = matrices[1];
	    
	    List<MotifPair> pairList = new ArrayList<MotifPair>();
	    for (int i = 0; i < set0.length; ++i) {
	    	for (int j = 0; j < set1.length; ++j) {
	    		double d = dMatrix.get(i, j);
	    		double f = fMatrix.get(i, j);
	    		double best;
	    		boolean bestIsFlipped;
	    		if (d < f) {
	    			best = d; bestIsFlipped = false;
	    		} else {
	    			best = f; bestIsFlipped = true;
	    		}
	    		if (best < threshold) {
	    			pairList.add(new MotifPair(set0[i], set1[j], best, bestIsFlipped));
	    		}
	    	}
	    }
	    return (MotifPair[]) pairList.toArray(new MotifPair[0]);
	}
    
    public MotifPair[] bestHits(Motif[] set0, Motif[] set1)
        throws Exception {
        Matrix2D[] matrices = getComparisonMatrix(set0, set1);
        Matrix2D dMatrix = matrices[0];
        Matrix2D fMatrix = matrices[1];
        
        List<MotifPair> pairList = bestHits(dMatrix, fMatrix, set0, set1, false);
                                            
        Collections.sort(pairList);
        Collections.reverse(pairList);
        return (MotifPair[]) pairList.toArray(new MotifPair[0]);
    }
    
    public MotifPair bestHit(Motif[] set0, Motif[] set1) 
    	throws Exception {
    	MotifPair[] pairs = bestHits(set0,set1);
    	
    	double minScore = Double.MAX_VALUE;
    	MotifPair bestPair = null;
    	for (MotifPair p : pairs) {
    		if (p.score < minScore) {
    			bestPair = p;
    		}
    	}
    	
    	return bestPair;
    }
    
    public MotifPair[] bestReciprocalHits(Motif[] set0, Motif[] set1)
        throws Exception {
        Matrix2D[] matrices = getComparisonMatrix(set0, set1);
        Matrix2D dMatrix = matrices[0];
        Matrix2D fMatrix = matrices[1];
        
        List<MotifPair> pairList = bestHits(dMatrix, fMatrix, set0, set1, false);
        pairList.retainAll(bestHits(new Transpose(dMatrix), new Transpose(fMatrix), set1, set0, true));
                                            
        Collections.sort(pairList);
        Collections.reverse(pairList);
        return (MotifPair[]) pairList.toArray(new MotifPair[0]);
    }
    
    private static class Transpose implements Matrix2D {
        private final Matrix2D raw;
        
        public Transpose(Matrix2D raw) {
            this.raw = raw;
        }
        
        public int rows() {
            return raw.columns();
        }
        public int columns() {
            return raw.rows();
        }
        
        public double get(int row, int column) {
            return raw.get(column, row);
        }
        
        public void set(int row, int column, double x) {
            raw.set(column, row, x);
        }
        
        public double[] getRaw() {
            throw new UnsupportedOperationException();
        }
    }

    protected static double div(Distribution d0, Distribution d1) 
	throws Exception
    {
        double cScore = 0.0;
            for (Iterator i = ((FiniteAlphabet) d0.getAlphabet()).iterator(); i.hasNext(); ) {
                 Symbol s= (Symbol) i.next();
                 cScore += Math.pow(d0.getWeight(s) - d1.getWeight(s), 2.0);
            }
        // return Math.sqrt(cScore);
        // return cScore;
        return Math.pow(cScore, 2.5 / 2.0);
    }
    
    public double compareMotifs(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1)
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

    
    public ScoreOffsetPair compareMotifsWithOffset
    	(WeightMatrix wm0, Distribution pad0, WeightMatrix wm1, Distribution pad1) 
    	throws Exception {
		double bestScore = Double.POSITIVE_INFINITY;
		int minPos = -wm1.columns();
		int maxPos = wm0.columns() + wm1.columns();
		int bestOffset = 0;
		
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
			if (score < bestScore) {
				bestScore = score;
				bestOffset =  offset;
			}
		}
		ScoreOffsetPair swf = new ScoreOffsetPair(bestOffset, bestScore);
		return swf;
	}
    
    public class ScoreOffsetPair implements Comparable<ScoreOffsetPair> {
    	private final int offset;
    	private final double score;
		
    	public ScoreOffsetPair(int offset, double score) {
    		this.offset = offset;
    		this.score = score;
    	}
    	
    	public int getOffset() {
			return offset;
		}
		public double getScore() {
			return score;
		}

		public int compareTo(ScoreOffsetPair o) {
			return Double.compare(this.offset, o.offset);
		}
    }
}
