package net.derkholm.nmica.extra.app.motif;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import net.derkholm.nmica.apps.MetaMotifFinder;
import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.maths.MathsTools;
import net.derkholm.nmica.model.Datum;
import net.derkholm.nmica.model.MultiplexContributionSampler;
import net.derkholm.nmica.model.metamotif.Dirichlet;
import net.derkholm.nmica.model.metamotif.MetaMotif;
import net.derkholm.nmica.model.metamotif.MetaMotifClippedSimplexPrior;
import net.derkholm.nmica.model.metamotif.MetaMotifColumnCloneSampler;
import net.derkholm.nmica.model.metamotif.MetaMotifIOTools;
import net.derkholm.nmica.model.metamotif.MetaMotifPrior;
import net.derkholm.nmica.model.metamotif.MetaMotifTools;
import net.derkholm.nmica.model.metamotif.bg.MetaMotifBackground;
import net.derkholm.nmica.model.metamotif.bg.MetaMotifDirichletBackground;
import net.derkholm.nmica.model.metamotif.bg.MetaMotifMixtureBackground;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifAlphaSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifMeanSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifPrecisionSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifRetrimSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifSeedSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifSlideSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifSymbolScalingSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifSymbolSwapSampler;
import net.derkholm.nmica.model.metamotif.sampler.MetaMotifZapSampler;
import net.derkholm.nmica.model.metamotif.sampler.SymbolWeightAlteringSampler;
import net.derkholm.nmica.model.motif.NMWeightMatrix;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import cern.colt.list.IntArrayList;

@App(overview="The NestedMICA meta motif sampler / spiking simulator.", generateStub=true)
@NMExtraApp(launchName="nmmetasim", vm=VirtualMachine.SERVER)
public class MetaMotifSimulator {
	//private File backgroundFile;
	private int num = -1;
	private int outputMotifLength = 20;
	private double[] freqs;
	private File[] inputMetaMotifFiles;
	private MetaMotif[] metaMotifs;
	private MetaMotifBackground background;
	private String outputFile;
	private String occOutputFile;
	
	private double[] backgroundAlphas;
	private int[] sampleNum;
	private int bgSampleNum = -1;
	private int[] maxNumMetaMotifHitsPerMotif;
	private double minAlphaSum = -1;
	private double maxAlphaSum = -1;
	private int randomMetaMotifs = -1;
	private double fixedAlphaSum = -1;
	private boolean randomizeColAlphaSums;
	private double[] fixedMeans = null;
	private double minMeanStep = -1;
	private double maxMeanStep = -1;
	private double minAlphaSumStep = -1;
	private double maxAlphaSumStep = -1;
	private int samplingChainLength = 1;
	private int minLength;
	private int maxLength;
	private int extraLength;
	private double minClip;
	private double maxClip;
	private boolean customSamplerSetup = true;
	private boolean useSlideSampler = false;
	private boolean useZapSampler = false;
	private boolean useMeanSampler = true;
	private boolean useAlphaSumSampler = true;
	private String[] useSamplerStrs;
	private MultiplexContributionSampler sampler;
	private File inputMotifDataFile;
	private double bgBiasAlpha;
	private double bgOtherAlpha;
	private double priorMeanSamplingAlphaSum;
	
	@Option(help="Background distribution alpha parameters (in the order A,C,G,T). Either this, -mixBG or -backgroundData has to be specified when spiking.", optional=true)
    public void setBackgroundAlphas(double[] ds) {
    	this.backgroundAlphas = ds;
    }
	
	@Option(help="Use a mixture background. Two arguments required: a 'biased' weight and the 'background' weight. " +
			"Either this, -backgroundAlphas or -backgroundData " +
			"has to be specified when spiking.", optional=true)
			
    public void setMixBG(double[] ds) {
		if (ds.length != 2) {
			System.err.println("When using -mixBG, specify the higher and the lower ");
			System.exit(1);
		} else if (ds[0] < 0 || ds[1] < 0) {
			System.err.println("-mixBG: positive values required for weights.");
			System.exit(1);
		}
    	this.bgBiasAlpha = ds[0];
    	this.bgOtherAlpha = ds[1];
    }
	
	@Option(help="Use the specified samplers on the input metamotifs (sample in metamotif space) : slide,zap,mean,alphasum " +
			"(note that mean and alphasum samplers " +
			"also require you to set bounds for the step size " +
			"with -alphaSumBounds and -meanSteBounds)", optional=true)
    public void setUseSamplers(String[] ss) {
		this.useSamplerStrs = ss;
    	for (String s : ss) {
    		if (s.equals("slide")) {
    			useSlideSampler = true;
    		}
    		else if (s.equals("zap")) {
    			useZapSampler = true;
    		}
    		else if (s.equals("mean")) {
    			useMeanSampler = true;
    		}
    		else if (s.equals("alphaSum")) {
    			useAlphaSumSampler = true;
    		}
    		else if (s.equals("all")) {
    			if (ss.length > 1) {
    				System.err.println("If -useSampler all specified, do not specify other samplers");
    				System.exit(1);
    			}
    			customSamplerSetup = false;
    		}
    	}
    	
    }
		
	@Option(help="Paths for XMS files containing input metamotifs", optional=true)
    public void setIn(File[] fs) {
    	this.inputMetaMotifFiles = fs;
    	
    }
	
	@Option(help="Path for an XMS files containing input motifs " +
			"(specify this option if you want to spike metamotifs " +
			"to an existing motif set file rather than to a background of Dirichlets, " +
			"i.e. this and -backgroundAlphas are mutually exclusive)", optional=true)
    public void setBackgroundData(File f) {
    	this.inputMotifDataFile = f;
    }
	
	@Option(help="Sample specified number(s) of motifs out of the metamotif(s). Number of counts has to match the number of metamotifs in the input metamotif file.", optional=true)
    public void setSampleNum(int[] is) {
    	this.sampleNum = is;
    }
	
	@Option(help="When sampling in metamotif space with -useSamplers, each sample taken initially is followed by a chain of further sampling operations (0 by default)", optional=true)
    public void setSampleChainLength(int i) {
    	this.samplingChainLength = i;
    }
	
	@Option(help="Sample specified number(s) of blank metamotifs out of the specified background distribution.", optional=true)
    public void setBgSampleNum(int i) {
    	this.bgSampleNum = i;
    }
	
	@Option(help="Output filename", optional=true)
    public void setOut(String s) {
    	this.outputFile = s;
    }
	
	@Option(help="Output filename for the occupancy information (when spiking metamotifs to simulated datasets)", optional=true)
    public void setOccOut(String s) {
    	this.occOutputFile = s;
    }
	
	@Option(help="Length of motif datasets to output (when spiking with -backgroundAlphas)", optional=true)
    public void setLength(int l) {
    	this.outputMotifLength = l;
    }
	
	@Option(help="Number of motifs to output", optional=true)
    public void setNum(int i) {
    	this.num  = i;
    }
	
	@Option(help="Sample random metamotifs " +
			"(used optionally in conjunction with " +
			"-randomColAlphaSums, -alphaSumBounds, -alphaSum, -fixedMeans", optional=true)
    public void setRandomMetaMotifs(int i) {
    	this.randomMetaMotifs = i;
    }
	
	@Option(help="Fix the generated random metamotif symbol frequencies " +
			"to the specified values " +
			"(the symbol for which each of the mean values is given is randomised)", optional=true)
    public void setFixedMeans(double[] ds) {
    	//TODO: Check validity of the input weight list
    	this.fixedMeans = ds;
    }
	
	@Option(help="Random alpha sum params for each column", optional=true)
    public void setRandomColAlphaSums(boolean b) {
    	this.randomizeColAlphaSums = b;
    }
	
	@Option(help="Minimum and maximum alpha sum bound", optional=true)
    public void setAlphaSumBounds(double[] ds) {
    	if (ds.length != 2) {
    		System.out.println("-alphaSumBounds: Incorrect number of arguments given: 2 required");
    		System.exit(1);
    	}
    	
    	this.minAlphaSum = ds[0];
    	this.maxAlphaSum = ds[1];
    }
	
	@Option(help="Mean scaling step bounds (two arguments needed: the min and the max)", optional=true)
    public void setMeanStepBounds(double[] ds) {
    	if (ds.length != 2) {
    		System.out.println("-meanStepBounds: Incorrect number of arguments given: 2 required");
    		System.exit(1);
    	}
    	
    	this.minMeanStep = ds[0];
    	this.maxMeanStep = ds[1];
    }
	
	@Option(help="Alpha sum scaling step bounds (two arguments needed: the min and the max)", optional=true)
    public void setAlphaSumStepBounds(double[] ds) {
    	if (ds.length != 2) {
    		System.out.println("-alphaSumScaleBounds: Incorrect number of arguments given: 2 required");
    		System.exit(1);
    	}
    	
    	this.minAlphaSumStep = ds[0];
    	this.maxAlphaSumStep = ds[1];
    }
	
	@Option(help="Minimum and maximum length for metamotif length distribution (two arguments needed: the min and the max)", optional=true)
    public void setLengthBounds(int[] is) {
    	if (is.length != 2) {
    		System.out.println("-lengthBounds: Incorrect number of arguments given: 2 required");
    		System.exit(1);
    	}
    	
    	this.minLength = is[0];
    	this.maxLength = is[1];
    }
	
	
	@Option(help="", optional=true)
    public void setAlphaSum(double d) {
    	this.fixedAlphaSum = d;
    }
	
	@Option(help="Relative frequencies of each of the input metamotif " +
			"(list of the same length " +
			"as there are metamotifs in the input file(s)", optional=true)
    public void setFreqs(double[] ds) {
    	this.freqs  = ds;
    }
	
	@Option(help="Maximum number of hits of metamotif(s) per data entry (motif). This should be a list of the same length as there are metamotifs.", optional=true)
    public void setMaxHits(int[] i) {
		for (int j : i)
		if (j <= 0) {
			System.err.println("Each value specified with -maxHits has to have a value of > 0");
			System.exit(1);
		}
    	this.maxNumMetaMotifHitsPerMotif = i;
    }
	
	public void main(String[] args) throws Exception {
		if (outputFile == null) {
			System.err.println("-out missing: you need to specify the output filename");
			System.exit(1);
		}
		
    	if (inputMetaMotifFiles != null)
    		metaMotifs = MetaMotifIOTools.loadMetaMotifsFromMultipleFiles(inputMetaMotifFiles);
		
    	/*
    	if (backgroundFile != null) {
    		MetaMotif[] backgroundMMs = MetaMotifIOTools
    									.loadMetaMotifSetXML(
    											new FileInputStream(backgroundFile));
    	
	    	if (backgroundMMs.length > 1) {
	    		System.err.println("More than one metamotif present in the input datafile");
	    		System.exit(1);
	    	}
	    	
	    	background = backgroundMMs[0].getColumn(0);
	    	backgroundMMs = null; //don't need this, let's GC it
	    	
    	} */
    	
    	if (backgroundAlphas != null) {
    		//TODO: Support other alphabets but DNA here
    		Dirichlet bDist = new MetaMotifDirichletBackground(1.0, DNATools.getDNA());
    		bDist.setWeight(DNATools.a(),backgroundAlphas[0]);
    		bDist.setWeight(DNATools.c(),backgroundAlphas[1]);
    		bDist.setWeight(DNATools.g(),backgroundAlphas[2]);
    		bDist.setWeight(DNATools.t(),backgroundAlphas[3]);
    		background = (MetaMotifDirichletBackground)bDist;
    	} else if (bgBiasAlpha > 0 && bgOtherAlpha > 0) {
    		MetaMotifMixtureBackground mixBg = new MetaMotifMixtureBackground(DNATools.getDNA(), bgBiasAlpha, bgOtherAlpha);
    		background = mixBg;
    	}
    	
    	if ((background == null) && (sampleNum == null) && this.inputMotifDataFile == null) {
    		if ((randomMetaMotifs <= 0) && (useSamplerStrs == null)) {
	    		System.err.println("Either -backgroundData, -backgroundAlphas, or -mixBG has to be specified.");
	    		System.exit(1);
    		}
    	} else if (sampleNum != null) {
    		if (background != null || inputMotifDataFile != null) {
    			System.err.println("Either -backgroundData|-backgroundAlphas|-mixBG, or " +
    						"-sampleNum has to be specified.");
    			System.exit(1);
    		}
    		
    		if (sampleNum.length != metaMotifs.length) {
    			System.err.println("Incorrect number of sampling counts specified.");
    			System.exit(1);
    		}
    		
    		sample();
    		System.exit(0);
    	}
    	
    	if (bgSampleNum < 0 && num > 0 && (background != null || inputMotifDataFile != null)) {
    		spike();
    	} else if (bgSampleNum > 0 && num < 0 && background != null) {
    		sampleBackground();
    	}
    	
    	if (this.randomMetaMotifs > 0) {
    		
    		makeRandomMetaMotifs();
    	}
    	
    	if (this.useSamplerStrs != null) {
    		
    		if (num <= 0) {
    			System.err.println("-num missing: you need to specify the number of samples drawn");
    			System.exit(1);
    		}
    		
    		MetaMotifIOTools.writeMetaMotifSetToMotifSetWithAnnotations(
    							new FileOutputStream(new File(outputFile)), 
    							sampleAroundMetaMotifs(metaMotifs));
    	}
    }

	private MetaMotif[] sampleAroundMetaMotifs(MetaMotif[] metaMotifs) throws Exception {
		if (sampler == null) setUpMultiplexSampler();
		int totalNum = metaMotifs.length * this.num;
		MetaMotif[] newMMs = new MetaMotif[totalNum];
		
		
		for (MetaMotif mm : metaMotifs)
			for (int i = 0; i < num; i++) {
				newMMs[i] = sampleAroundMetaMotif(mm);
				mm.setName(mm.getName() + "_sample" + i);
				System.out.print(".");
			}
		System.out.println("\n");
		return newMMs;
	}

	private void sampleBackground() throws Exception {
		Motif[] motifs = new Motif[bgSampleNum];
		for (int c = 0; c < bgSampleNum; c++) {
			Motif m = new Motif();
			m.setName("bg_sample_" + c);
			WeightMatrix wm = sampleWeightMatrixFromDirichlet(background, outputMotifLength);
			m.setWeightMatrix(wm);
			m.setThreshold(0);
			motifs[c] = m;
		}
		
		MotifIOTools.writeMotifSetXML(new FileOutputStream(outputFile), motifs);
	}
	
	public class MotifSetModificationMaskMatrix {
		private int[][] modificationMask;
		private HashMap<MetaMotif, Integer> metaMotifIds;
		private int[][] metaMotifCountsPerTargetMotif;
		private MetaMotif[] metaMotifs;
				
		private void setUpMetaMotifIds() {
			for (int i = 0; i < metaMotifs.length; i++)
				metaMotifIds.put(metaMotifs[i], i);
		}
		
		public MotifSetModificationMaskMatrix(int numMotifs, Motif[] outputMotifs, MetaMotif[] metaMotifs) {
			modificationMask = new int[numMotifs][];
			for (int i = 0; i < modificationMask.length; i++) {
				int cols = outputMotifs[i].getWeightMatrix().columns();
				modificationMask[i] = new int[cols];
				for (int j = 0; j < modificationMask[i].length; j++)
					modificationMask[i][j] = -1;
			}
				
			this.metaMotifs  = metaMotifs;
			metaMotifCountsPerTargetMotif = new int[numMotifs][metaMotifs.length];
			metaMotifIds = new HashMap<MetaMotif, Integer>();
			
			setUpMetaMotifIds();
		}
		
		//TODO: This isn't working
		public boolean modificationPossible(int row, int startOffset, int hitLength) {
			int[] modificationMaskRow = modificationMask[row];
			for (int i = startOffset; i < (startOffset + hitLength); i++)
				if (modificationMaskRow[i] >= 0) return false; //a modification already done at this position
			return true;
		}
		
		public void addModification(MetaMotif mm, int row, int startOffset, int hitLength) {
			int[] modificationMaskRow = modificationMask[row];
			int metaMotifId = metaMotifIds.get(mm);
			
			//TODO: Check that this + 1 here is correct
			for (int i = startOffset; i < (startOffset + hitLength + 1); i++) {
				if (i < 0) continue;
				if (i >= modificationMaskRow.length) continue;
				
				modificationMaskRow[i] = metaMotifId;
			}
			
			metaMotifCountsPerTargetMotif[row][metaMotifIds.get(mm)]++;
			
		}
		
		public int longestUnmodifiedStretchLength(int row) {
			int[] modificationMaskRow = modificationMask[row];
			int max = 0;
			
			int curStretch = 0;
			for (int i = 0, length = modificationMaskRow.length; i < length; i++)
				if (modificationMaskRow[i] == -1)
					curStretch++;
				else {
					if (max < curStretch) max = curStretch;
					curStretch = 0;
				}
			
			if (max < curStretch)
				max = curStretch;
				
			return max;
		}
		
		public int overallLongestUnmodifiedStretchLength() {
			int max = 0;
			
			for (int i = 0; i < modificationMask.length; i++) {
				int thisMax = longestUnmodifiedStretchLength(i);
				if (max < thisMax)
					max = thisMax;
			}
			return max;
		}
		
		public int numHits(int row, MetaMotif mm) {
			return metaMotifCountsPerTargetMotif[row][metaMotifIds.get(mm)];
		}
		
		
		//TODO: Support an XML format too 
		//(either an XML snippet or extra annotations could be put in an XMS file)
		public String toOccupancyStr() {
			StringBuffer strBuf = new StringBuffer();
			for (int m = 0; m < modificationMask.length; m++) {
				IntArrayList indices = new IntArrayList();
				for (int mm = 0; mm < modificationMask[m].length; mm++) {
					//-1 is the value for an empty location
					if (!indices.contains(modificationMask[m][mm]) &! (modificationMask[m][mm] == -1))
						indices.add(modificationMask[m][mm]);
				}
				
				for (int mm = 0; mm < indices.size(); mm++) {
					strBuf.append(indices.get(mm));
					if (mm < (indices.size() - 1))
						strBuf.append(ModelScoreEvaluator.OCC_COL_SEPARATOR);
				}
				strBuf.append("\n");
			}	
			return strBuf.toString();
		}
	}
	
	private List<Symbol> symbols(FiniteAlphabet alphabet) {
		List<Symbol> symbols = new ArrayList<Symbol>(alphabet.size());
		
		for (Iterator it = alphabet.iterator(); it.hasNext();)
			symbols.add((Symbol) it.next());
		
		return symbols;
	}
	
	private void setUpMultiplexSampler() {
		MetaMotifSymbolScalingSampler meanSampler = null;
		MetaMotifSymbolScalingSampler alphaSumSampler = null;
		MetaMotifSlideSampler slideSampler = null;
		MetaMotifZapSampler zapSampler = null;
		
		MetaMotifClippedSimplexPrior mmPrior = new MetaMotifClippedSimplexPrior(
				DNATools.getDNA(), 
				minLength, 
				maxLength, 
				maxLength + extraLength, 
				minClip, 
				maxClip,
				this.minAlphaSum,
				this.maxAlphaSum,
				priorMeanSamplingAlphaSum);

		if (this.customSamplerSetup) {
			MultiplexContributionSampler mcsSampler = new MultiplexContributionSampler();
			if (this.useMeanSampler) {
				if (minMeanStep > 0) {
					meanSampler = new MetaMotifSymbolScalingSampler(
							minMeanStep, 
							maxMeanStep, 
							MetaMotifSymbolScalingSampler.ExpectedValueSamplingPolicy.MULTIPLY);
				} else { 
					System.err.println("When specifying sampler \"mean\", you must also specify -meanStepBounds");
					System.exit(1);
				}
			}
			if (this.useAlphaSumSampler) {
				if (maxAlphaSumStep > 0) {
					alphaSumSampler = new MetaMotifSymbolScalingSampler(maxAlphaSumStep, 
							MetaMotifSymbolScalingSampler.ExpectedValueSamplingPolicy.MULTIPLY);
				} else {
					System.err.println("When specifying sampler \"alphasum\", you must also specify -alphaSumStepBounds");
					System.exit(1);
				}
			}

			if (this.useSlideSampler) {
				slideSampler = new MetaMotifSlideSampler();
			}

			if (this.useZapSampler) {
				zapSampler = new MetaMotifZapSampler(mmPrior);
			}

			if (meanSampler != null)
				mcsSampler.addSampler(meanSampler, 1.0);
			if (alphaSumSampler != null)
				mcsSampler.addSampler(alphaSumSampler, 1.0);
			if (slideSampler != null)
				mcsSampler.addSampler(slideSampler, 1.0);
			if (zapSampler != null)
				mcsSampler.addSampler(zapSampler, 1.0);
			
			sampler = mcsSampler;
		} else {
		//set up in the same way as in MetaMotifFinder
		sampler = MetaMotifSimulator.setUpMultiplexContributionSampler(
				null, 
				mmPrior, 
				this.maxMeanStep, 
				this.maxAlphaSumStep, 
				-1, //maxSeedVariance, not used as we're not seeding from existing data
				minLength, maxLength);
		}
	}
	
	private MetaMotif sampleAroundMetaMotif(MetaMotif mm) throws Exception {
		return (MetaMotif)sampler.sample(mm, null).getVariate();
	}
	
	private void makeRandomMetaMotifs() throws FileNotFoundException, Exception {
		MetaMotif[] metaMotifs = new MetaMotif[this.randomMetaMotifs];
		
		for (int i = 0; i < this.randomMetaMotifs; i++) {
		
			MetaMotif mm;
			
			if (this.fixedAlphaSum > 0)
				mm = MetaMotifSimulator.sampleRandomMetaMotifWithPrecision(this.outputMotifLength, fixedAlphaSum);
				
			else if ((0 < minAlphaSum) && (minAlphaSum < maxAlphaSum))
				if (this.randomizeColAlphaSums)
					mm = MetaMotifSimulator.sampleRandomMetaMotifWithRandomPrecisionsForEachColumn(this.outputMotifLength, minAlphaSum, maxAlphaSum);
				else 
					mm = MetaMotifSimulator.sampleRandomMetaMotifWithRandomPrecision(this.outputMotifLength, minAlphaSum, maxAlphaSum);
			else {
				System.err.println("Need to either supply a positive value for -alphaSum or two valid positive values -alphaSumBounds");
				System.exit(1);
				return;
			}
			
			//TODO: This is a bit of a hackery (weight modified for an unnecessary cycle)
			if (this.fixedMeans != null) {
				for (int j = 0; j < mm.columns(); j++) {
					List<Symbol> sampledSymbols = symbols(mm.getAlphabet());
					double alphaSum = mm.getColumn(j).alphaSum();
					int symIndex = 0;
					
					do {
						int index = MathsTools.randomInt(sampledSymbols.size());
						Symbol sym = sampledSymbols.remove(index);
						mm.getColumn(j).setWeight(sym, fixedMeans[symIndex++] * alphaSum);
					} while (sampledSymbols.size() > 0);
				}
			}
			
			
			metaMotifs[i] = mm;
		}
		
		MetaMotifIOTools.writeMetaMotifSetToMotifSetWithAnnotations(new FileOutputStream(new File(this.outputFile)), metaMotifs);
		
	}

	private void spike() throws IllegalSymbolException,
			IllegalAlphabetException, Exception, FileNotFoundException {

    	if (outputFile == null) {
    		System.exit(1);
    		System.err.println("-out missing (output filename)");
    	}
    	
    	if (freqs == null) {
    		System.err.println("-freqs missing");
    		System.exit(1);
    	} else {
    		for (double f : freqs)
    			if (f < 0) {
    				System.err.println("Invalid frequency given: " + f);
    				System.exit(1);
    			}
    	}
    	
    	//set the default values in case nothing was specified (i.e. no limitations on hit numbers)
    	if (maxNumMetaMotifHitsPerMotif == null) {
    		maxNumMetaMotifHitsPerMotif = new int[metaMotifs.length];
    		for (int i = 0; i < maxNumMetaMotifHitsPerMotif.length; i++) {
    			maxNumMetaMotifHitsPerMotif[i] = Integer.MAX_VALUE;
    		}
    	} else if (maxNumMetaMotifHitsPerMotif.length != metaMotifs.length) {
    		System.err.println("The correct number of values required for -maxHits is equal to the number of metamotifs specified (" + metaMotifs.length + ")");
    		System.exit(1);
    		
    	}

		Motif[] outputMotifs = null;
    	if (inputMotifDataFile != null && background != null) {
    		System.err.println("Specify either -backgroundAlphas or -backgroundData when spiking");
    		System.exit(1);
    	} else if (inputMotifDataFile != null) {
    		outputMotifs = MotifIOTools.loadMotifSetXML(new FileReader(inputMotifDataFile));
    	} else {
    		//simulated background
    		outputMotifs = new Motif[num];
    		for (int i = 0; i < outputMotifs.length; i++) {
        		Motif m = new Motif();
        		m.setWeightMatrix(sampleWeightMatrixFromDirichlet(background, outputMotifLength));
        		m.setName("motif" + i);
        		//TODO: Also output the distribution parameters as an annotation to the file
        		outputMotifs[i] = m;
        	}
    	}
    	
		MotifSetModificationMaskMatrix modificationMask = 
			new MotifSetModificationMaskMatrix(num, outputMotifs, metaMotifs);
		
    	for (int i = 0; i < metaMotifs.length; i++) {
    		MetaMotif mm = metaMotifs[i];
    		double freq = freqs[i];
    		
    		//TODO: You should check and fail if it's not realistic to fit the specified count in there
    		int count = (int)Math.round(freq  * num);
    		
    		//these are randomised within the do--while loop
    		int targetIndex;
    		int offset;
    		WeightMatrix targetWM;
    		
    		for (int c = 0; c < count; c++) {
    			WeightMatrix wm = MetaMotifTools.sampleFromMetaMotif(mm);
    			
    			int attempts = 0;
    			boolean maxNumHitsMet = true;
    			do {
	    			//TODO: Support only placing a max number of metamotif smaples per motif
					targetIndex = MathsTools.randomInt(outputMotifs.length);
					targetWM = outputMotifs[targetIndex].getWeightMatrix();
	    			int maxOffset = targetWM.columns() - wm.columns() - 1;
    			
	    			if (maxOffset < 0)
    						throw new IllegalArgumentException(
    						"Simulated output motif length is too short " +
    						"to fit samples of the weight matrix");
    			
    				offset = MathsTools.randomInt(maxOffset);
	    			//TODO: Output the metamotif locations as annotations 
	    			//(can be then visualised with mxplor)
    				attempts++;
    				
    				//TODO: You should cache results of overallLongestUnmodifiedStretchLength() inside the modification mask
    				//TODO: Will run indefinitely if ovreallLongestUnmodifiedStretchLongth i
    				//TODO: Treat == separately (no randomisation needed then)
    				maxNumHitsMet = modificationMask.numHits(targetIndex, mm) >= maxNumMetaMotifHitsPerMotif[i];
    				if (modificationMask.longestUnmodifiedStretchLength(targetIndex) >= mm.columns()) {
    					continue;
    				} else if (modificationMask.overallLongestUnmodifiedStretchLength() >= mm.columns()) {
    					continue;
    				} else {
    					System.err.println("Could not spike with " + mm.getName() + " (too small number / too small length of simulated background motifs");
    					System.exit(1);
    				}
    					
    			} while (!modificationMask.modificationPossible(targetIndex,offset,mm.columns()) || maxNumHitsMet);
    			
    			modifyWeights(targetWM, wm, offset);
    			modificationMask.addModification(mm, targetIndex, offset, mm.columns());
    		}
    	}
    	
    	MotifIOTools.writeMotifSetXML(new FileOutputStream(outputFile), outputMotifs);
    	
    	if (occOutputFile != null) {
    		BufferedWriter writer = 
    			new BufferedWriter(
    					new OutputStreamWriter(
    							new FileOutputStream(new File(occOutputFile))));
    		
    		String occStr = modificationMask.toOccupancyStr();
    		System.out.println(occStr);
    		writer.write(occStr);
    		writer.close();
    	}
	}

	private void sample() throws Exception, FileNotFoundException {
		int totalCount = 0;
		for (int i : sampleNum) totalCount += i;
		
		Motif[] outputMotifs = new Motif[totalCount];
		int index = 0;
		for (int c = 0; c < sampleNum.length; c++) {
			MetaMotif mm = metaMotifs[c];
			int count = sampleNum[c];
			for (int i = 0; i < count; i++) {
				String name = mm.getName() + "_sample"+i;
				outputMotifs[index++] = MetaMotifTools.sampleMotifFromMetaMotif(mm, name);
			}
		}
		
		MotifIOTools.writeMotifSetXML(new FileOutputStream(outputFile), outputMotifs);
	}

	private static MetaMotif sampleRandomMetaMotifWithPrecision(int length, double alphaSum) {
		Dirichlet[] dists = new Dirichlet[length];
		
		for (int i = 0; i < length; i++)
			dists[i] = MetaMotif.sampleUniformly(alphaSum);
		
		MetaMotif mm;
		try {
			mm = new MetaMotif(dists);
		} catch (IllegalSymbolException e) {
			e.printStackTrace();
			throw new BioError("Assertion failed : illegal symbol exception caught when making a random metamotif");
		}
		
		return mm;
	}
	
	private static MetaMotif sampleRandomMetaMotifWithRandomPrecisionsForEachColumn(int length, double minAlphaSum, double maxAlphaSum) {
		MetaMotif mm = sampleRandomMetaMotifWithPrecision(length,1);
		
		for (int i = 0 ; i < mm.columns(); i++) {
			Dirichlet dist = mm.getColumn(i);
			
			for (Iterator it = ((FiniteAlphabet) dist.getAlphabet())
					.iterator(); it.hasNext();) {
				Symbol sym = (Symbol) it.next();
				try {
					dist.setWeight(sym, dist.getWeight(sym) * (minAlphaSum + (maxAlphaSum - minAlphaSum) * Math.random()));
				} catch (IllegalSymbolException e) {
					throw new BioError("Illegal symbol exception caught", e);
				}
			}
		}
		
		return mm;
	}
	
	
	private static MetaMotif sampleRandomMetaMotifWithRandomPrecision(int length, double minAlphaSum, double maxAlphaSum) {
		return sampleRandomMetaMotifWithPrecision(length, minAlphaSum + (maxAlphaSum - minAlphaSum) * Math.random());
	}
	
	private void modifyWeights(WeightMatrix targetWM, WeightMatrix wm, int startOffset) throws IllegalSymbolException {
		for (int j = 0; j < wm.columns(); j++) {
			if ((startOffset + j) < 0) continue;
			if ((startOffset + j) >= wm.columns()) continue;
			
			Distribution d = targetWM.getColumn(startOffset + j);
			Distribution newD = wm.getColumn(j);
			for (Iterator it = ((FiniteAlphabet) d.getAlphabet())
					.iterator(); it.hasNext();) {
				Symbol sym = (Symbol) it.next();
				d.setWeight(sym,newD.getWeight(sym));
			}
		}
	}

	private WeightMatrix sampleWeightMatrixFromDirichlet
	(MetaMotifBackground distribution, int metaMotifLength) {
		Distribution[] dists = new Distribution[metaMotifLength];
		for (int j = 0; j < metaMotifLength; j++)
			dists[j] = distribution.sampleDistribution();

		try {
			return new NMWeightMatrix(dists,dists.length,0);
		} catch (IllegalAlphabetException e) {
			e.printStackTrace();
			throw new BioError("Assertion failed : illegal alphabet caught", e);
		}
	}
	
	private WeightMatrix sampleWeightMatrixFromDirichlet
			(Dirichlet distribution, int metaMotifLength) {
		Distribution[] dists = new Distribution[metaMotifLength];
		for (int j = 0; j < metaMotifLength; j++)
			dists[j] = distribution.sampleDistribution();
		
		try {
			return new NMWeightMatrix(dists,dists.length,0);
		} catch (IllegalAlphabetException e) {
			e.printStackTrace();
			throw new BioError("Assertion failed : illegal alphabet caught", e);
		}
	}

	// TODO: this is only used in SpikeMetaMotifs, move it there
	public static MultiplexContributionSampler setUpMultiplexContributionSampler(
			Datum[] data, MetaMotifPrior mmPrior, double alphaMaxScaleStep,
			double alphaSumMaxScale, double maxSeedVariance, int minLength,
			int maxLength) {
		MultiplexContributionSampler mcsSampler = new MultiplexContributionSampler();
	
		/*
		mcsSampler
				.addSampler(
						new MetaMotifSymbolScalingSampler(
								alphaMaxScaleStep,
								MetaMotifSymbolScalingSampler.ExpectedValueSamplingPolicy.MULTIPLY),
						MEAN_SCALING_SAMPLER_WEIGHT);
		*/
		
		mcsSampler.addSampler(new MetaMotifMeanSampler(alphaMaxScaleStep),MetaMotifFinder.getMeanScalerWeight());
		
		mcsSampler.addSampler(new MetaMotifAlphaSampler(alphaMaxScaleStep, SymbolWeightAlteringSampler.SamplingPolicy.GAUSSIAN), MetaMotifFinder.getAlphaScalerWeight());
		
		/*
		mcsSampler
				.addSampler(
						new MetaMotifSymbolScalingSampler(
								alphaSumMaxScale,
								MetaMotifSymbolScalingSampler.PrecisionSamplingPolicy.MULTIPLY),
						ALPHA_SUM_SCALING_SAMPLER_WEIGHT);*/
		
		mcsSampler.addSampler(new MetaMotifPrecisionSampler(alphaSumMaxScale, SymbolWeightAlteringSampler.SamplingPolicy.GAUSSIAN), MetaMotifFinder.getAlphaSumScalerWeight());
	
		/*
		mcsSampler
				.addSampler(new MetaMotifSymbolScalingSampler(MetaMotifSymbolScalingSampler.SymbolWeightSwappingPolicy.SWAP_RANDOMLY),
						SYMBOL_SWAP_SAMPLER_WEIGHT);
		*/
		
		mcsSampler.addSampler(new MetaMotifSymbolSwapSampler(DNATools.getDNA(), 
				MetaMotifSymbolSwapSampler.SymbolWeightSwappingPolicy.SWAP_RANDOMLY),
				MetaMotifFinder.getSymbolSwapSamplerWeight());
		
		mcsSampler
				.addSampler(new MetaMotifSlideSampler(), MetaMotifFinder.getSlideSamplerWeight());
	
		if (mmPrior != null) {
			mcsSampler.addSampler(new MetaMotifZapSampler(mmPrior),
					MetaMotifFinder.getZapSamplerWeight());
			
			if (minLength != maxLength) {
				mcsSampler.addSampler(new MetaMotifRetrimSampler(mmPrior),
					MetaMotifFinder.getRetrimSamplerWeight());
			}
		}
	
		if (data != null)
			mcsSampler
					.addSampler(
							new MetaMotifSeedSampler(
									mmPrior,
									data,
									maxSeedVariance,
									MetaMotifSeedSampler.SeedPrecisionPolicy.UNIFORM_RANDOM_SEED_PRECISION,
									MetaMotifSeedSampler.SeedLengthPolicy.MAX_SEED_LENGTH),
							MetaMotifFinder.getSeedSamplerWeight());
	
		mcsSampler.addSampler(new MetaMotifColumnCloneSampler(),
				MetaMotifFinder.getColumnCloneSamplerWeight());
	
		return mcsSampler;
	}
}
