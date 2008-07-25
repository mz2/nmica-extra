package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import biobits.utils.Collects;
import biobits.utils.Function;
import biobits.utils.Function2;

import cern.jet.stat.Gamma;


@App(overview="Compare real and theoretical motif-score distributions", generateStub=true)
@NMExtraApp(launchName="nmhistcomp", vm=VirtualMachine.SERVER)
public class MotifMatchHistogramComparitor {
	private double bucketSize = 1.0;
	private double confidence = 0.01;
	private File theoreticalFile;
	private File realFile;
	
	@Option(help="Bucket size in bits (default=1.0)", optional=true)
	public void setBucketSize(double bucketSize) {
		this.bucketSize = bucketSize;
	}

	@Option(help="Confidence threshold (default=0.01)", optional=true)
	public void setConfidence(double confidence) {
		this.confidence = confidence;
	}
	
	@Option(help="Theoretical distribution (output format of nmweightwords)")
	public void setTheoretical(File f) {
		this.theoreticalFile = f;
	}
	
	@Option(help="Observed distribution in input sequences " +
			"as tab separated records with sequence identifier and score per line " +
			"(e.g. first and sixth column of a GFF)")
	public void setObserved(File f) {
		this.realFile = f;
	}

	private List<HEl> bucket(File f)
		throws Exception {
		List<HEl> bl = new ArrayList<HEl>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			StringTokenizer t = new StringTokenizer(line);
			t.nextToken();
			int bucket = (int) Math.floor(Math.abs(Double.parseDouble(t.nextToken())) / bucketSize);
			double weight = 1.0;
			if (t.hasMoreTokens()) {
				weight = Double.parseDouble(t.nextToken());
			}
			bl.add(new HEl(bucket, weight));
		}
		return bl;
	}
	
	private static final Function2<Integer, Integer, Integer> MAX = new Function2<Integer, Integer, Integer>() {
		public Integer apply(Integer param1, Integer param2) {
			return Math.max(param1, param2);
		}
	};
	
	private Integer max(List<HEl> l) {
		Function<HEl, Integer> b = new Function<HEl, Integer>() {
			public Integer apply(HEl param1) {
				return param1.bucket;
			}
		};
		return Collects.reduce(Collects.map(b, l), MAX);
	}
	
	private double[] hist(Iterable<HEl> i, int buckets) {
		double[] h = new double[buckets];
		for (HEl el : i) {
			h[el.bucket] += el.weight;
		}
		return h;
	}
	
	public void main(String[] args) throws Exception {
		List<HEl> reference = bucket(theoreticalFile);
		List<HEl> real = bucket(realFile);
		int buckets = 1 + Math.max(max(reference), max(real));
		double[] refHist = hist(reference, buckets);
		double refTot = 0;
		for (double r : refHist) {
			refTot += r;
		}
		double[] realHist = hist(real, buckets);
		for (int b = 0; b < buckets; ++b) {
			double refFrac = refHist[b] / refTot;
			double mode = (1.0 * realHist[b]) / real.size();
			
			int alpha = (int) realHist[b] + 1;
            int beta = real.size() - (int) realHist[b] + 1;
            
            double barMin = 0;
            double barMax = 1.0;
            while (barMin < mode && Gamma.incompleteBeta(alpha, beta, barMin) < confidence) {
                barMin += 0.00001;
            }
            while (barMax > mode && Gamma.incompleteBeta(alpha, beta, barMax)> (1.0 - confidence)) {
                barMax -= 0.00001;
            }
            
            System.out.printf("%d\t%g\t%g\t%g\t%g%n", b, refFrac, mode, barMin, barMax);
		}
	}
	
	private static class HEl {
		public final int bucket;
		public final double weight;
		
		HEl(int bucket, double weight) {
			this.bucket = bucket;
			this.weight = weight;
		}
	}
}
