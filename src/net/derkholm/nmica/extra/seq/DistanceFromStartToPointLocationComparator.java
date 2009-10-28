package net.derkholm.nmica.extra.seq;
import java.util.Comparator;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.PointLocation;


public class DistanceFromStartToPointLocationComparator implements Comparator<StrandedFeature> {
	private int referencePoint;
	
	public DistanceFromStartToPointLocationComparator(PointLocation refPos) {
		this.referencePoint = refPos.getMin();
	}

	private int distance(StrandedFeature s, int referencePoint) {
		int ref,refStart,refEnd;
		if (s.getStrand().equals(StrandedFeature.POSITIVE)) {
			ref = Math.abs(s.getLocation().getMin() - referencePoint);
		} else if (s.getStrand().equals(StrandedFeature.NEGATIVE)) {
			ref = Math.abs(s.getLocation().getMax() - referencePoint);
		} else {
			refStart = Math.abs(s.getLocation().getMin() - referencePoint);
			refEnd = Math.abs(s.getLocation().getMax() - referencePoint);
			ref = Math.min(refStart, refEnd);
		}
		return ref;
	}
	
	public int compare(StrandedFeature s0, StrandedFeature s1) {
		
		int s0Dist = distance(s0,this.referencePoint);
		int s1Dist = distance(s1,this.referencePoint);
		
		if (s0Dist < s1Dist) return -1;
		else if (s0Dist > s1Dist) return 1;
		
		return 0;
	}
}