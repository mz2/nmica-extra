package net.derkholm.nmica.motif.align;

import net.derkholm.nmica.motif.Motif;

public class MotifPairWithOffset extends MotifPair {
    protected final int offset;
    
	public MotifPairWithOffset(Motif m1, Motif m2, double score, boolean flip,int offset) {
		super(m1, m2, score, flip);
		this.offset = offset;
	}

}
