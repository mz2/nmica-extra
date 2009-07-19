package net.derkholm.nmica.extra.app;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
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

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.maths.MathsTools;
import net.derkholm.nmica.model.motif.extra.ScoredHit;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import biobits.utils.IOTools;

@App(overview="Area under an ROC", generateStub=true)
@NMExtraApp(launchName = "nmrocauc")
public class MotifHitROCAUCalculator {
	private int bootstraps = 10000;
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
	private TreeSet<String> motifNames;
	
	@Option(help="Permute labels again (test)", optional=true)
	public void setPermuteLabels(boolean b) {
		this.test = b;
	}
	
	@Option(help=
			"Input score list is in E-value format, " +
			"but use raw score anyway", 
			optional=true)
	public void setEvalsRaw(boolean b) {
		this.evalsRaw = b;
	}
	
	@Option(help="Input score list is in E-value format", 
			optional=true)
	public void setEvals(boolean b) {
		this.evals = b;
	}
	
	@Option(help="Target motif whose hits to seek (hits to other motifs are ignored)", optional=true)
	public void setTarget(String s) {
		this.target = s;
	}
	
	@Option(help="The number of bootstraps to make (default=10000)", optional=true)
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
	
	public void setPositiveHits(List<ScoredHit> hits) {
		this.positiveHits = filter(hits, true);
	}
	
	public void setNegativeHits(List<ScoredHit> hits) {
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
		throws Exception
	{
		List<ScoredHit> positiveHits = read(new File(args[0]), true);
		List<ScoredHit> negativeHits = read(new File(args[1]), false);
		
		this.motifNames = new TreeSet<String>();
		for (ScoredHit hit : positiveHits) {
			motifNames.add(hit.getMotifName());
		}
		for (ScoredHit hit : negativeHits) {
			motifNames.add(hit.getMotifName());
		}
		
		motifPositiveHitMap = mapHitsToMotifs(positiveHits, motifNames);
		motifNegativeHitMap = mapHitsToMotifs(negativeHits, motifNames);
		
		List<MotifROCAUCSummary> summaries = new ArrayList<MotifROCAUCSummary>();
		
		for (String motifName : motifNames) {
			List<ScoredHit> l = new ArrayList<ScoredHit>();
			l.addAll(motifPositiveHitMap.get(motifName));
			l.addAll(motifNegativeHitMap.get(motifName));
			
			int numTrue = 0, numFalse = 0;
			for (ScoredHit s : l) {
				if (s.isPositive()) {
					++numTrue;
				} else {
					++numFalse;
				}
			}
			
			Collections.shuffle(l); // destabilize;
			Collections.sort(l, new Comparator<ScoredHit>() {
				public int compare(ScoredHit arg0, ScoredHit arg1) {
					return MathsTools.sign(arg1.getScore() - arg0.getScore());
				}
			});
			
			if (test) {
				Collections.shuffle(l);
			}
			double auc = rocAuc2(l, numTrue, numFalse);
			int over = 0;
			for (int b = 0; b < bootstraps; ++b) {
				Collections.shuffle(l);
				if (rocAuc2(l, numTrue, numFalse) >= auc) {
					++over;
				}
			}
			summaries.add(new MotifROCAUCSummary(motifName, auc, (1.0 * over) / bootstraps));
		}
		
		Collections.sort(summaries);
		for (MotifROCAUCSummary sum : summaries) {
			System.out.printf(
				"%s\t%g\t%g%n", 
				sum.getMotifName(),
				sum.getAuc(),
				sum.getBootstrapFraction());
		}
	}
	
	private double rocAuc(Collection<ScoredHit> l, int tot) {
		int numSeen = 0;
		int numTrue = 0;
		double auc = 0.0;
		for (ScoredHit sh : l) {
			++numSeen;
			if (sh.isPositive()) {
				++numTrue;
				auc += (1.0 / tot) * ((1.0 * numTrue) / numSeen);
			}
		}
		return  auc;
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
	
	public class MotifROCAUCSummary implements Comparable {
		private String motifName;
		private double auc;
		private double bootstrapFraction;

		public MotifROCAUCSummary(String motifName, double auc, double bootstrapFraction) {
			this.motifName = motifName;
			this.auc = auc;
			this.bootstrapFraction = bootstrapFraction;
		}

		public String getMotifName() {
			return motifName;
		}

		public void setMotifName(String motifName) {
			this.motifName = motifName;
		}

		public double getAuc() {
			return auc;
		}

		public void setAuc(double auc) {
			this.auc = auc;
		}

		public double getBootstrapFraction() {
			return bootstrapFraction;
		}

		public void setBootstrapFraction(double bootstrapFraction) {
			this.bootstrapFraction = bootstrapFraction;
		}

		public int compareTo(Object o) {
			if (!(o instanceof MotifROCAUCSummary)) {
				throw new IllegalArgumentException(
					"Trying to compare a MotifROCAUCSummary instance against a " + 
					o.getClass().getCanonicalName());
			}
			
			MotifROCAUCSummary sum = (MotifROCAUCSummary) o;
			return Double.compare(this.getAuc(), sum.getAuc());
		}
	}
}
