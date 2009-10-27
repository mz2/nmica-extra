package net.derkholm.nmica.extra.peak;

public class Peak {
	private String refSeqName;
	private int min;
	private int max;

	public Peak(String refSeq, int i) {
		this.refSeqName = refSeq;
		this.min = i;
	}

	public String getRefSeqName() {
		return refSeqName;
	}

	public void setRefSeqName(String refSeqName) {
		this.refSeqName = refSeqName;
	}

	public int getMin() {
		return min;
	}

	public void setMin(int min) {
		this.min = min;
	}

	public int getMax() {
		return max;
	}

	public void setMax(int max) {
		this.max = max;
	}
}