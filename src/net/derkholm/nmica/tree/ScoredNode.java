package net.derkholm.nmica.tree;

public class ScoredNode extends Node {
	protected double score;
	public ScoredNode(double score) {
		super();
		this.score = score;
	}

	public ScoredNode(Object data, double score) {
		super(data);
		this.score = score;
	}


	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
}
