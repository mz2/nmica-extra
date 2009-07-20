package net.derkholm.nmica.model.analysis;

public class BucketComparisonElement {

	private int bucket;
	private double referenceFraction;
	private double mode;
	private double barMin;
	private double barMax;

	public BucketComparisonElement(int bucket, double referenceFraction,
			double mode, double barMin, double barMax) {
		this.bucket = bucket;
		this.referenceFraction = referenceFraction;
		this.mode = mode;
		this.barMin = barMin;
		this.barMax = barMax;
	}

	public int getBucket() {
		return bucket;
	}

	public void setBucket(int bucket) {
		this.bucket = bucket;
	}

	public double getReferenceFraction() {
		return referenceFraction;
	}

	public void setReferenceFraction(double referenceFraction) {
		this.referenceFraction = referenceFraction;
	}

	public double getMode() {
		return mode;
	}

	public void setMode(double mode) {
		this.mode = mode;
	}

	public double getBarMin() {
		return barMin;
	}

	public void setBarMin(double barMin) {
		this.barMin = barMin;
	}

	public double getBarMax() {
		return barMax;
	}

	public void setBarMax(double barMax) {
		this.barMax = barMax;
	}
}