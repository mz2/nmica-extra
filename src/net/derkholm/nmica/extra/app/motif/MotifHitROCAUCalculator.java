package net.derkholm.nmica.extra.app.motif;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.maths.MathsTools;
import net.derkholm.nmica.model.analysis.MotifROCAUCSummary;
import net.derkholm.nmica.model.analysis.ScoredHit;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import biobits.utils.IOTools;

@App(overview="Area under an ROC", generateStub=true)
@NMExtraApp(launchName = "nmrocauc")
public class MotifHitROCAUCalculator {
	private int bootstraps = 50000;
	private String target = null;
	private Set<String> whiteList = null;
	private Set<String> blackList = null;
	private boolean evals = false;
	private boolean evalsRaw = false;
	private boolean test = false;
	private List<ScoredHit> positiveHits;
	private List<ScoredHit> negativeHits;
	private HashMap<String,List<ScoredHit>> motifPositiveHitMap = new HashMap<String,List<ScoredHit>>();
	private HashMap<String,List<ScoredHit>> motifNegativeHitMap = new HashMap<String,List<ScoredHit>>();
	private Collection<String> motifNames;
	private File motifs;
	private File positiveSeqs;
	private File negativeSeqs;
	private File positives;
	private File negatives;
	private File outputFile;
	
	@Option(help="Permute labels again (test)", optional=true)
	public void setPermuteLabels(boolean b) {
		this.test = b;
	}
	
	@Option(help="Input score list is in E-value format: motif \\t seq \\t eval", 
			optional=true)
	public void setEvalsRaw(boolean b) {
		this.evalsRaw = b;
	}
	
	@Option(help="Input score list is in E-value format: motif \\t seq \\t maxscore \\t eval", 
			optional=true)
	public void setEvals(boolean b) {
		this.evals = b;
	}
	
	@Option(help="Target motif whose hits to seek (hits to other motifs are ignored)", optional=true)
	public void setTarget(String s) {
		this.target = s;
	}
	
	@Option(help="The number of bootstraps to make (default=50000)", optional=true)
	public void setBootstraps(int b) {
		this.bootstraps = b;
	}
	
	@Option(help="White list of sequence identifiers (allowed hits in the positive sequence set)", optional=true)
	public void setWhiteList(Reader r) 
		throws IOException
	{
		whiteList = new HashSet<String>();
		BufferedReader br = new BufferedReader(r);
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			whiteList.add(line);
		}
	}
	
	@Option(help="Black list of sequence identifiers (allowed hits in the negative sequence set)", optional=true)
	public void setBlackList(Reader r) 
		throws IOException
	{
		blackList = new HashSet<String>();
		BufferedReader br = new BufferedReader(r);
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			blackList.add(line);
		}
	}
	
	@Option(help="Motif set file (optional - you can either give input with -motifs and -seqs or with )", optional=true)
	public void setMotifs(File f) 
		throws Exception {
		this.motifs = f;
	}

	@Option(help="Positive sequences to score", optional=true)
	public void setPositiveSeqs(File seqs) {
		this.positiveSeqs = seqs;
	}

	@Option(help="Negative sequences to score", optional=true)
	public void setNegativeSeqs(File seqs) {
		this.negativeSeqs = seqs;
	}

	@Option(help="Positive hits. Default input format: motif \\t seq \\t maxscore \\t eval", optional=true)
	public void setPositive(File f) {
		this.positives = f;
	}
	
	@Option(help="Negative hits. Default input format: motif \\t seq \\t maxscore \\t eval", optional=true)
	public void setNegative(File f) {
		this.negatives = f;
	}
	
	@Option(help="Output file")
	public void setOut(File f) {
		this.outputFile = f;
	}
	
	private int threads = 1;
	private ExecutorService threadPool;

	@Option(help="Number of threads (default = 1)", optional = true)
	public void setThreads(int threads) {
		if (threads < 1) {
			System.err.println("-threads needs to be >= 1");
			System.exit(1);
		}
		this.threads = threads;
	}
	
	private void setPositiveHits(List<ScoredHit> hits) {
		this.positiveHits = filter(hits, true);
	}

	private void setNegativeHits(List<ScoredHit> hits) {
		this.negativeHits = filter(hits, false);
	}
	
	
	
	private List<ScoredHit> filter(List<ScoredHit> hits, boolean label) {
		List<ScoredHit> filteredHits = new ArrayList<ScoredHit>();
		
		for (ScoredHit hit : hits) {
			String seq = hit.getSeqName();
			String motif = hit.getMotifName();
			double score = hit.getScore();
			if (evals) {
				double eval = hit.getEValue();
				if (eval <= 0) {
					score = 100000;
				} else {
					score = -Math.log(eval) / Math.log(10);
				}
			}
			
			boolean accepted = false;
			if (whiteList == null) {
				accepted = true;
			} else if (whiteList.contains(seq)) {
				accepted = true;
			} else {
				int uso = seq.indexOf('_');
				accepted = (uso > 0 && whiteList.contains(seq.substring(0, uso)));
			}
			if (blackList != null && blackList.contains(seq)) {
				accepted = false;
			}
			if (accepted) {
				if (target == null || target.equals(motif)) {
					filteredHits.add(new ScoredHit(motif, seq, label, score, score));
				}
			}
		}
		
		return filteredHits;
	}
	
	private List<ScoredHit> read(File f, boolean label)
		throws Exception {
		
		List<ScoredHit> hits = new ArrayList<ScoredHit>();
		
		BufferedReader br = IOTools.fileBufferedReader(f);
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			StringTokenizer toke = new StringTokenizer(line);
			String seq = toke.nextToken();
			String motif = toke.nextToken();
			double score = Double.parseDouble(toke.nextToken());
			if (evals) {
				{
					String tmp = motif;
					motif = seq;
					seq = tmp;
				}
				
				double eval = Double.parseDouble(toke.nextToken());
				if (eval <= 0) {
					score = 100000;
				} else {
					score = -Math.log(eval) / Math.log(10);
				}
			} else if (evalsRaw) {
				{
					String tmp = motif;
					motif = seq;
					seq = tmp;
				}
			}
			
			boolean accepted = false;
			if (whiteList == null) {
				accepted = true;
			} else if (whiteList.contains(seq)) {
				accepted = true;
			} else {
				int uso = seq.indexOf('_');
				accepted = (uso > 0 && whiteList.contains(seq.substring(0, uso)));
			}
			if (blackList != null && blackList.contains(seq)) {
				accepted = false;
			}
			if (accepted) {
				if (target == null || target.equals(motif)) {
					hits.add(new ScoredHit(motif, seq, label, score, score));
				}
			}
		}
		return hits;
	}
	
	public HashMap<String,List<ScoredHit>> mapHitsToMotifs(List<ScoredHit> hits, Collection<String> motifNames) {
		HashMap<String,List<ScoredHit>> map = new HashMap<String, List<ScoredHit>>();
		
		for (String motifName : motifNames) {
			map.put(motifName, new ArrayList<ScoredHit>());
			for (ScoredHit hit : hits) {
				if (hit.getMotifName().equals(motifName)) {
					map.get(motifName).add(hit);
				}
			}			
		}
		return map;
	}
	
	public void main(String[] args)
		throws Exception {
		
		PrintStream os = null;
		if (this.outputFile == null) {
			os = System.out;
		} else {
			os = new PrintStream(new FileOutputStream(this.outputFile, true));
		}
		
		if ((motifs != null) && (positiveSeqs != null) && (negativeSeqs != null)) {
			System.err.println("Calculating e-values...");
			MotifSetEmpiricalEValueCalculator eValueCalc = new MotifSetEmpiricalEValueCalculator();
			eValueCalc.setThreads(threads);
			eValueCalc.setBootstraps(this.bootstraps);
			eValueCalc.setCollectHits(true);
			
			eValueCalc.setMotifs(new FileReader(motifs));
			
			eValueCalc.setPositiveHits(true);
			eValueCalc.setSeqs(positiveSeqs);
			System.err.println("Calculating positive hits...");
			eValueCalc.calculate();
			System.err.println("Filtering positive hits...");
			this.setPositiveHits(new ArrayList<ScoredHit>(eValueCalc.collectedHits()));
			
			eValueCalc.clearCollectedHits();
			
			eValueCalc.setPositiveHits(false);
			eValueCalc.setSeqs(negativeSeqs);
			System.err.println("Calculating negative hits...");
			eValueCalc.calculate();
			System.err.println("Filtering negative hits...");
			this.setNegativeHits(new ArrayList<ScoredHit>(eValueCalc.collectedHits()));
			
		} else if (positives != null && negatives != null) {
			System.err.println("Filtering positive hits...");
			this.setPositiveHits(read(positives, true));
			System.err.println("Filtering negative hits...");
			this.setNegativeHits(read(negatives, false));			
		} else {
			System.err.println(
				"You have to either specify -motifs,-positiveSeqs,-negativeSeqs OR -positives,-negatives");
			System.exit(1);
		}
		

		if (motifs == null) {
			this.motifNames = new TreeSet<String>();
			for (ScoredHit hit : positiveHits) {
				motifNames.add(hit.getMotifName());
			}
			for (ScoredHit hit : negativeHits) {
				motifNames.add(hit.getMotifName());
			}
		} else {
			List<String> motifNs = new ArrayList<String>();
			Motif[] motifSet = MotifIOTools.loadMotifSetXML(new FileReader(motifs));
			for (Motif m : motifSet) motifNs.add(m.getName());
			motifNames = motifNs;
		}
		int motifCount = motifNames.size();
		
		motifPositiveHitMap = mapHitsToMotifs(positiveHits, motifNames);
		motifNegativeHitMap = mapHitsToMotifs(negativeHits, motifNames);
		
		List<Future<MotifROCAUCSummary>> summaryFutures = new ArrayList<Future<MotifROCAUCSummary>>();
		
		if (threads > motifCount) {
			System.err.println("Specified number of threads is larger than the number of motifs in the input. " +
					"Will use only " + motifCount + " threads.");
		}
		
		threadPool = Executors.newFixedThreadPool(threads);
		
		for (String motifName : motifNames) {
			List<ScoredHit> l = new ArrayList<ScoredHit>();
			l.addAll(motifPositiveHitMap.get(motifName));
			l.addAll(motifNegativeHitMap.get(motifName));
			
			summaryFutures.add(threadPool.submit(new ROCAUCTask(motifName, l, bootstraps, test)));
		}
		
		List<MotifROCAUCSummary> summaries = new ArrayList<MotifROCAUCSummary>();
		for (Future<MotifROCAUCSummary> sum : summaryFutures) {summaries.add(sum.get());}
		Collections.sort(summaries);
		
		for (MotifROCAUCSummary sum : summaries) {
			os.printf(
				"%s\t%g\t%g%n", 
				sum.getMotifName(),
				sum.getAuc(),
				sum.getBootstrapFraction());
		}
		
		threadPool.shutdown();
	}
	
	private static class ROCAUCTask implements Callable<MotifROCAUCSummary> {
		private String motifName;
		private int bootstraps;
		private final List<ScoredHit> hits;
		private boolean test;
		
		public ROCAUCTask(String motifName, List<ScoredHit> hits, int numBootstraps, boolean test) {
			this.hits = hits;
			this.motifName = motifName;
			this.bootstraps = numBootstraps;
			this.test = test;
		}
		
		
		public MotifROCAUCSummary call() throws Exception {
			
			
			int numTrue = 0, numFalse = 0;
			for (ScoredHit s : hits) {
				if (s.isPositive()) {
					++numTrue;
				} else {
					++numFalse;
				}
			}
			
			Collections.shuffle(hits); // destabilize;
			Collections.sort(hits, new Comparator<ScoredHit>() {
				public int compare(ScoredHit arg0, ScoredHit arg1) {
					return MathsTools.sign(arg1.getScore() - arg0.getScore());
				}
			});
			
			if (test) {
				Collections.shuffle(hits);
			}
			
			double auc = rocAuc2(hits, numTrue, numFalse);
			int over = 0;
			for (int b = 0; b < bootstraps; ++b) {
				Collections.shuffle(hits);
				if (rocAuc2(hits, numTrue, numFalse) >= auc) {
					++over;
				}
			}
			return new MotifROCAUCSummary(motifName, auc, (1.0 * over) / bootstraps);
		}
		
		private double rocAuc2(Collection<ScoredHit> l, int trues, int falses) {
			int numTrue = 0;
			int numFalse = 0;
			double auc = 0.0;
			for (ScoredHit sh : l) {
				if (sh.isPositive()) {
					++numTrue;
					auc += (1.0 / trues) * ((1.0 * numFalse) / falses);
				} else {
					++numFalse;
				}
			}
			return  1.0 - auc;
		}
	}
	
	
}
