package net.derkholm.nmica.extra.peak;

public class Window {
	protected String refSeqName;
	protected int beginCoord;
	protected int endCoord;
	protected double depth;
	protected double pvalue;
	protected double gcContent;
	
	public Window(
			String refSeqName, 
			int beginCoord, 
			int endCoord, 
			double depth, double pvalue) {
		this.refSeqName = refSeqName;
		this.beginCoord = beginCoord;
		this.endCoord = endCoord;
		this.depth = depth;
		this.pvalue = pvalue;
	}

	public String getRefSeqName() {
		return refSeqName;
	}

	public void setRefSeqName(String refSeqName) {
		this.refSeqName = refSeqName;
	}

	public int getBeginCoord() {
		return beginCoord;
	}

	public void setBeginCoord(int beginCoord) {
		this.beginCoord = beginCoord;
	}

	public int getEndCoord() {
		return endCoord;
	}

	public void setEndCoord(int endCoord) {
		this.endCoord = endCoord;
	}

	public double getDepth() {
		return depth;
	}

	public void setDepth(double depth) {
		this.depth = depth;
	}

	public double getPvalue() {
		return pvalue;
	}

	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}

	public double getGcContent() {
		return gcContent;
	}

	public void setGcContent(double gcContent) {
		this.gcContent = gcContent;
	}
}