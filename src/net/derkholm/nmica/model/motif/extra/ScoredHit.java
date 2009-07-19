package net.derkholm.nmica.model.motif.extra;

public class ScoredHit {
	private final String motifName;
	private final String seqName;
	private final boolean positive;
	private final double score;
	private final double eValue;
	
	public ScoredHit(String motifName, 
						String seqName, 
						boolean label, 
						double score,
						double eValue) {
		this.motifName = motifName;
		this.seqName = seqName;
		this.positive = label;
		this.score = score;
		this.eValue = eValue;
	}

	public String getMotifName() {
		return motifName;
	}

	public String getSeqName() {
		return seqName;
	}

	public boolean isPositive() {
		return positive;
	}

	public double getScore() {
		return score;
	}

	public double getEValue() {
		return eValue;
	}
}