package net.derkholm.nmica.extra.seq.nextgen;

import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class SAMPileup {
	private int extendedLength;
	private short[] pileup;

	public SAMPileup(String refName, int refLength, int extendedLength) {
		this.pileup = new short[refLength];
		this.extendedLength = extendedLength;
	}
	
	public void add(SAMRecord rec) {
		int start = rec.getAlignmentStart();
		int end = rec.getAlignmentEnd();
		if (!rec.getReadNegativeStrandFlag()) {for (int i = start; i < end; i++) {pileup[i]++;}}
		else {for (int i = end - extendedLength; i < end; i++) {pileup[i]++;}}
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