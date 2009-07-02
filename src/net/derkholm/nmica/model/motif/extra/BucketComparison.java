package net.derkholm.nmica.model.motif.extra;

import java.util.List;

import net.derkholm.nmica.model.HistogramElementIFace;

import biobits.utils.Collects;
import biobits.utils.Function;
import biobits.utils.Function2;
import cern.jet.stat.Gamma;

public class BucketComparison {

	protected double[] referenceHistogram;
	protected List<HistogramElementIFace> referenceHits;
	protected double[] observedHistogram;
	protected List<HistogramElementIFace> observedHits;
	
	//lazily fetched
	protected double referenceTotal = Double.NEGATIVE_INFINITY;
	private int numBuckets;
	private double bucketSize;
	
	public BucketComparison(
			List<HistogramElementIFace> referenceHits, 
			List<HistogramElementIFace> observedHits,
			double bucketSize) {
		this.numBuckets = 1 + Math.max(max(referenceHits), max(observedHits));
		this.bucketSize = bucketSize;
		this.setReferenceHits(referenceHits);
		this.setObservedHits(observedHits);
		
		//System.err.println("Number of buckets: " + this.numBuckets);
	}
	
	public double referenceTotal() {
		if (referenceTotal == Double.NEGATIVE_INFINITY) {
			double refTot = 0;
			for (double r : referenceHistogram) {
				refTot += r;
			}
			referenceTotal = refTot;
		}
		return referenceTotal;
	}
	
	public int alpha(int bucket) {
		return (int) observedHistogram[bucket] + 1;
	}
	
	public int beta(int bucket) {
		return observedHits.size() - (int) observedHistogram[bucket] + 1;
	}

	public double mode(int bucket) {
		return (1.0 * observedHistogram[bucket]) / observedHits.size();
	}
	

	public double barMin(int bucket, double confidence) {
		double barMin = 0;
		//System.err.println("mode:" + mode(bucket));
		//System.err.println("alpha:" + alpha(bucket));
		//System.err.println("beta:" + beta(bucket));
		
		while (barMin < mode(bucket) && 
				(Gamma.incompleteBeta(
				this.alpha(bucket), this.beta(bucket), barMin) < confidence)) {
			barMin += 0.00001;
		}
		return barMin;
	}

	public double barMax(int bucket, double confidence) {
		double barMax = 1.0;

		//System.err.println("bucket: " + bucket + " confidence: " + confidence);
		while (barMax > this.mode(bucket) && 
				(Gamma.incompleteBeta(alpha(bucket), 
									 beta(bucket), barMax) > (1.0 - confidence))) {
			barMax -= 0.00001;
		}
		return barMax;
	}
	
	public double referenceFraction(int bucket) {
		return referenceHistogram[bucket] / this.referenceTotal();
	}

	public String outputString(int bucket, double confidence) {
		
		
		return String.format("%d\t%g\t%g\t%g\t%g%n", 
				bucket, 
				referenceFraction(bucket), 
				mode(bucket), 
				barMin(bucket, confidence), 
				barMax(bucket, confidence));
	}
	
	public BucketComparisonElement compare(int bucket, double confidence) {
		return new BucketComparisonElement(
					bucket, 
					referenceFraction(bucket), 
					mode(bucket), 
					barMin(bucket, confidence), 
					barMax(bucket, confidence));
	}

	public List<HistogramElementIFace> getReferenceHits() {
		return referenceHits;
	}

	private void setReferenceHits(List<HistogramElementIFace> referenceHits) {
		this.referenceHits = referenceHits;
		this.referenceHistogram = hist(this.referenceHits, numBuckets);
		//System.err.println("Reference histogram:");
		//for (double d : this.referenceHistogram) {
		//	System.err.println(d);
		//}
	}

	public List<HistogramElementIFace> getObservedHits() {
		return observedHits;
	}

	private void setObservedHits(List<HistogramElementIFace> observedHits) {
		this.observedHits = observedHits;
		this.observedHistogram = hist(this.observedHits, numBuckets);
		//System.err.println("Observed histogram:");
		//for (double d : this.observedHistogram) {
		//	System.err.println(d);
		//}
	}
	
	private static final Function2<Integer, Integer, Integer> MAX = new Function2<Integer, Integer, Integer>() {
		public Integer apply(Integer param1, Integer param2) {
			return Math.max(param1, param2);
		}
	};
	
	private Integer max(List<HistogramElementIFace> l) {
		if (l.size() == 0) {
			return new Integer(0);
		}
		Function<HistogramElementIFace, Integer> b = new Function<HistogramElementIFace, Integer>() {
			public Integer apply(HistogramElementIFace param1) {
				return param1.bucket();
			}
		};
		return Collects.reduce(Collects.map(b, l), MAX);
	}
	
	private double[] hist(Iterable<HistogramElementIFace> i, int buckets) {
		double[] h = new double[buckets];
		for (HistogramElementIFace el : i) {
			h[el.bucket()] += el.hitWeight();
		}
		return h;
	}
	
	public int buckets() {
		return numBuckets;
	}

}
