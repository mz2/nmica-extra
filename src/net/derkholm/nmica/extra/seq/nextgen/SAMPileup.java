package net.derkholm.nmica.extra.seq.nextgen;

import net.sf.samtools.SAMRecord;

import org.biojava.bio.BioError;

public class SAMPileup {
	private int extendedLength;
	private short[] pileup;
	private int refLength;
	private String refName;

	public SAMPileup(String refName, int refLength, int extendedLength) {
		this.pileup = new short[refLength];
		this.extendedLength = extendedLength;
		this.refLength = refLength;
		this.refName = refName;
	}

	public void add(SAMRecord rec) {
		if (!rec.getReferenceName().equals(this.refName)) {
			throw new BioError(String.format("Unexpected reference sequence %s (expecting %s)",rec.getReferenceName(), this.refName));
		}

		int start = rec.getAlignmentStart();
		int end = rec.getAlignmentEnd();
		if (!rec.getReadNegativeStrandFlag()) {
			for (int i = Math.max(start,0), endpos = Math.min(this.refLength, start + extendedLength); i < endpos; i++) {pileup[i]++;}
		} else {
			for (int i = Math.max(0, Math.min(this.refLength,end - extendedLength)); i < end; i++) {pileup[i]++;}
		}
	}

	public int depthAt(int i) {
		return pileup[i];
	}

	public int getExtendedLength() {
		return extendedLength;
	}

	public void setExtendedLength(int extendedLength) {
		this.extendedLength = extendedLength;
	}
}