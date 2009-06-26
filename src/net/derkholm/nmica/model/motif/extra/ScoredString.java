package net.derkholm.nmica.model.motif.extra;

import net.derkholm.nmica.extra.app.ScoredSequenceHit;

public class ScoredString implements ScoredSequenceHit {
	
	protected String string;
	protected double score;
	protected double bgScore;
	
	public ScoredString(String string, double score) {
		this.string = string;
		this.score = score;
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
		this.bgScore = bgScore;
	}

	public double score() {
		return score;
	}

	public double weight() {
		return bgScore;
	}
	
	
}
