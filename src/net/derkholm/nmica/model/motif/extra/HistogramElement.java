package net.derkholm.nmica.model.motif.extra;

import net.derkholm.nmica.model.HistogramElementIFace;

public class HistogramElement implements HistogramElementIFace {
	private int bucket;
	private double weight;
	
	public HistogramElement(int bucket, double weight) {
		this.bucket = bucket;
		this.weight = weight;
	}
	
	public int bucket() {
		return bucket;
	}
	
	public double hitWeight() {
		return weight;
	}
}