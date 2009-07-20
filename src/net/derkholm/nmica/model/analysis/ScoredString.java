package net.derkholm.nmica.model.analysis;

import net.derkholm.nmica.model.ScoredSequenceHit;


public class ScoredString implements ScoredSequenceHit {
	
	protected String string;
	protected double score;
	protected double bgScore;
	protected int bucket;
	
	public ScoredString(String string, double score) {
		this.string = string;
		this.score = Math.abs(score);
	}

	public String getString() {
		return string;
	}

	public void setString(String string) {
		this.string = string;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public double getBgScore() {
		return bgScore;
	}

	public void setBgScore(double bgScore) {
		this.bgScore = Math.abs(bgScore);
	}

	public double score() {
		return score;
	}

	public double hitWeight() {
		return Math.abs(bgScore);
	}

	public void setBucket(int bucket) {
		this.bucket = bucket;
	}

	public int bucket() {
		return this.bucket;
	}
	
	
}
