package net.derkholm.nmica.extra.app.motifs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.maths.NativeMath;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;
import net.derkholm.nmica.utils.CollectTools;

import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.DoubleAlphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SimpleAtomicSymbol;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.bio.symbol.DoubleAlphabet.DoubleSymbol;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Output background likelihoods of a sequence given a mosaic background model", generateStub=true)
@NMExtraApp(launchName="nmmostool", vm=VirtualMachine.SERVER)
public class MosaicTool {
	private File mosaic;
	
	private boolean logLikelihoods = true;
	private boolean viterbi = false;
	
	@Option(help="A NestedMICA mosaic background XML file")
	public void setMosaic(File mosaic) {
		this.mosaic = mosaic;
	}
	
	@Option(help="Output log likelihoods (default=true)", optional=true)
	public void setLog(boolean b) {
		this.logLikelihoods = b;
	}
	
	@Option(help="Output the most likely state path (default=false)", optional=true)
	public void setViterbi(boolean b) {
		this.viterbi = b;
	}
	
	private final static Pattern labelPattern = Pattern.compile("patch([0-9]+)");
	private final static Pattern statePattern = Pattern.compile("\\(([\\w,\\s]+)\\)");
	private final static Pattern logLikelihoodPattern = Pattern.compile("\\((-{0,1}\\d+\\.{0,1}\\d*)");
	
	
	
	public void main(String[] args)
		throws Exception {
		XMLInputFactory factory = XMLInputFactory.newInstance();
		XMLStreamReader r = factory.createXMLStreamReader(new FileInputStream(this.mosaic));
		Mosaic m = MosaicIO.readMosaic(r);
		MosaicSequenceBackground  mosaic = new MosaicSequenceBackground(m.getDistributions(), m.getTransition());
		
		for (String fn : args) {
			if (!(new File(fn).exists())) {
				System.err.println("File " + fn + " does not exist.");
				System.exit(1);
			}
		}
		
		Distribution[] background = mosaic.getBackgroundDistributions();
        int backgroundOrder = background[0].getAlphabet().getAlphabets().size();
		
		
		for (String fn : args) {
			SequenceIterator iter = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(fn)));
			
			while (iter.hasNext()) {
				Sequence seq = iter.nextSequence();
								
				if (!viterbi) {
					double[] likelihoods = bgLikelihood(mosaic, seq, logLikelihoods);

					//System.out.println("bglikelihood:" + likelihoods.length);
					//System.out.println("seqlength   :" + seq.length());
					
					for (int i = 1 + backgroundOrder; i < seq.length(); i++) {
						System.out.println(seq.subStr(i, i)+"\t"+likelihoods[i-1]);
					}
				} else {
					DP dp = mosaic.getBackgroundDP();
					
					StatePath sp = dp.viterbi(new SymbolList[]{
							SymbolListViews.orderNSymbolList(seq, backgroundOrder)}, ScoreType.PROBABILITY);
					
					List states = sp.toList();
					
					for (Object o : states) {
						AtomicSymbol sym = (AtomicSymbol)o;
						List syms = sym.getSymbols();
						if (syms.size() != 3) {
							throw new BioException("Expected three basis symbols but got " + syms.size());
						} else {
							DoubleAlphabet.DoubleSymbol likelihood = (DoubleSymbol) syms.get(0);
							SimpleEmissionState state = (SimpleEmissionState) syms.get(1);
							SimpleAtomicSymbol symbols = (SimpleAtomicSymbol) syms.get(2);
							
							List seqSyms = symbols.getSymbols();
							
							Object symO = seqSyms.get(backgroundOrder-1);
							Symbol bsym = (Symbol) symO;
							
							String symN = ((FiniteAlphabet)seq.getAlphabet()).getTokenization("token").tokenizeSymbol(bsym);
							System.out.println(symN + "\t" + likelihood.getName() + "\t" + state.getName());
						}
					}
				}
				
				System.out.println("-");
			}
		}
		
        /*
        int[] seqLens = new int[background.length];
        
        
        SymbolList seq = inseq;
        if (backgroundOrder > 1) {
            seq = SymbolListViews.orderNSymbolList(seq, backgroundOrder);
        } 
        

        Matrix2D post = new SimpleMatrix2D(1 + seq.length(), background.length);
        if (chunk <= 0) {
        	dpChunk(mosaic, post, seq, 1, seq.length(), 1);
        } else {
        	int globMinTake = 1;
        	while (globMinTake < seq.length()) {
        		int globMaxTake = Math.min(globMinTake + chunk - 1, seq.length());
        		int globMinDP = Math.max(1, globMinTake - chunkOverlap);
            	int globMaxDP = Math.min(globMaxTake + chunkOverlap, seq.length());
            	
            	System.err.printf("Chunking %d-%d%n", globMinDP, globMaxDP);
            	dpChunk(mosaic, post, seq.subList(globMinDP, globMaxDP), globMinTake - globMinDP + 1, globMaxTake - globMinDP + 1, globMinTake);
            	
            	globMinTake = globMaxTake + 1;
        	}
        	
        }*/
    }

	private double[] bgLikelihood(MosaicSequenceBackground mosaic, SymbolList seq, boolean log) 
		throws IllegalSymbolException, IllegalAlphabetException {
		double[] tempHood = new double[seq.length()];
        Location seqMask = mosaic.backgroundSymbolLikelihood(seq, tempHood);
        
        AlphabetIndex index = AlphabetManager
        						.getAlphabetIndex((FiniteAlphabet)seq.getAlphabet());
        
        List<Integer> indexList = new ArrayList<Integer>();
        List<Double> hoodList = new ArrayList<Double>();
        double bgHood = 0;
        
        for (int pos = 1; pos <= seq.length(); ++pos) {
            if (seqMask.contains(pos)) {
                hoodList.add(new Double(tempHood[pos - 1]));
                bgHood += tempHood[pos - 1];
                Symbol s = seq.symbolAt(pos);
                if (s instanceof AtomicSymbol) {
                    indexList.add(new Integer(index.indexForSymbol(s)));
                } else {
                    indexList.add(new Integer(-1));
                }
            }
        }
        
        double[] symbolHoodLog = CollectTools.toDoubleArray(hoodList);	        
        
        if (log) { 
        	return symbolHoodLog;
        } else {
	        double[] symbolHood = new double[symbolHoodLog.length];
	        for (int i = 0; i < symbolHood.length; ++i) {
	            symbolHood[i] = NativeMath.exp2(symbolHoodLog[i]);
	        }
	        return symbolHood;
        }
        
	}
	
	/*
	
	private void dpChunk(MosaicSequenceBackground mosaic, Matrix2D post, SymbolList seq, int minTake, int maxTake, int postOffset)
		throws Exception {
		DP dp = mosaic.getBackgroundDP();
		Distribution[] background = mosaic.getBackgroundDistributions();
        int backgroundOrder = background[0].getAlphabet().getAlphabets().size();
        
        SingleDPMatrix forwardMatrix = (SingleDPMatrix) dp.forwardMatrix(
            new SymbolList[] {seq},
            ScoreType.PROBABILITY
        );
        System.err.printf("Forward = %g%n", forwardMatrix.getScore());
        double score = forwardMatrix.getScore();
        SingleDPMatrix backwardMatrix = (SingleDPMatrix) dp.backwardMatrix(
            new SymbolList[] {seq},
            ScoreType.PROBABILITY
        );
        System.err.printf("Backward = %g%n", backwardMatrix.getScore());
        
        int prevWhichMax = -1;
        int whichMax = -1;
        int consecutiveHits = 0;
        SymbolList symList = new SimpleSymbolList(DNATools.getDNA());
        
        for (int pos = minTake; pos <= maxTake; ++pos) {
            Symbol sym = seq.symbolAt(pos);
            State[] states = forwardMatrix.states();
            double[] fcol = forwardMatrix.scores[pos];
            double[] bcol = backwardMatrix.scores[pos];
            double sh = 0;
            
            double max = Double.NEGATIVE_INFINITY;
            
            int position = pos - minTake + postOffset;
            

            prevWhichMax = whichMax;
            whichMax = -1;
            
            System.out.println(states.length);
            for (int s = 0; s < states.length; ++s) {
                Matcher nameMatcher = namePattern.matcher(states[s].getName());
                if (nameMatcher.matches()) {
                    int state = Integer.parseInt(nameMatcher.group(1));
                    double weight = Math.exp(fcol[s] + bcol[s] - score);
                    
                    if (weight > max) {
                    	max = weight;
                    	whichMax = state;
                    }
                }
            }
            
            System.out.println("Class:");
    		System.out.println(sym.getClass());
    		System.out.println(sym.getMatches().getClass());
    		System.out.println(sym.getName());
    		
            if (prevWhichMax == whichMax) {
            	consecutiveHits++;
        		symList.edit(new Edit(symList.length()+1,DNATools.getDNA(),sym));
            } else {
            	if (consecutiveHits >= hitThreshold) {
            		System.out.println(
        				DNATools.getDNA().getTokenization("token").tokenizeSymbolList(symList));
            		
            		symList = new SimpleSymbolList(DNATools.getDNA());
            	}
            	
            	consecutiveHits = 0;
            }
        }
	} */
}
