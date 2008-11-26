package net.derkholm.nmica.extra.app;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.maths.MathsTools;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import biobits.utils.IOTools;

@App(overview="Area under an ROC", generateStub=true)
@NMExtraApp(launchName = "nmrocauc")
public class MotifHitROCAUCalculator {
	private int bootstraps = 1000;
	private String target = null;
	private Set<String> whiteList = null;
	private Set<String> blackList = null;
	private boolean evals = false;
	private boolean evalsRaw = false;
	private boolean test = false;
	
	@Option(help="...", optional=true)
	public void setTest(boolean b) {
		this.test = b;
	}
	
	@Option(help="Input score list is in E-value format, but use raw score anyway", optional=true)
	public void setEvalsRaw(boolean b) {
		this.evalsRaw = b;
	}
	
	@Option(help="Input score list is in E-value format", optional=true)
	public void setEvals(boolean b) {
		this.evals = b;
	}
	
	@Option(help="...", optional=true)
	public void setTarget(String s) {
		this.target = s;
	}
	
	@Option(help="...", optional=true)
	public void setBootstraps(int b) {
		this.bootstraps = b;
	}
	
	@Option(help="...", optional=true)
	public void setWhiteList(Reader r) 
		throws IOException
	{
		whiteList = new HashSet<String>();
		BufferedReader br = new BufferedReader(r);
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			whiteList.add(line);
		}
	}
	
	@Option(help="...", optional=true)
	public void setBlackList(Reader r) 
		throws IOException
	{
		blackList = new HashSet<String>();
		BufferedReader br = new BufferedReader(r);
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			blackList.add(line);
		}
	}
	
	private void read(Collection<? super ScoredHit> l , File f, boolean label)
		throws Exception
	{
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
					l.add(new ScoredHit(label, score));
				}
			}
		}
	}
	
	public void main(String[] args)
		throws Exception
	{
		List<ScoredHit> l = new ArrayList<ScoredHit>();
		read(l, new File(args[0]), true);
		read(l, new File(args[1]), false);
		int numTrue = 0, numFalse = 0;
		for (ScoredHit s : l) {
			if (s.isTrue) {
				++numTrue;
			} else {
				++numFalse;
			}
		}
		
		Collections.shuffle(l); // destabilize;
		Collections.sort(l, new Comparator<ScoredHit>() {
			public int compare(ScoredHit arg0, ScoredHit arg1) {
				return MathsTools.sign(arg1.score - arg0.score);
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
		System.out.printf("%g\t%g%n", auc, (1.0 * over) / bootstraps);
	}
	
	private double rocAuc(Collection<ScoredHit> l, int tot) {
		int numSeen = 0;
		int numTrue = 0;
		double auc = 0.0;
		for (ScoredHit sh : l) {
			++numSeen;
			if (sh.isTrue) {
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
			if (sh.isTrue) {
				++numTrue;
				auc += (1.0 / trues) * ((1.0 * numFalse) / falses);
			} else {
				++numFalse;
			}
		}
		return  1.0 - auc;
	}
	
	private static class ScoredHit {
		public final boolean isTrue;
		public final double score;
		
		public ScoredHit(boolean b, double s) {
			this.isTrue = b;
			this.score = s;
		}
	}
}
