package net.derkholm.nmica.extra.seq;
import java.util.Comparator;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.PointLocation;


public class DistanceToPointLocationComparator implements Comparator<StrandedFeature> {
	private PointLocation referencePosition;

	public DistanceToPointLocationComparator(PointLocation refPos) {
		this.referencePosition = refPos;
	}

	public int compare(StrandedFeature o1, StrandedFeature o2) {
		int ref0Start, ref0End;
		ref0Start = Math.abs(o1.getLocation().getMin() - referencePosition.getMin());
		ref0End = Math.abs(o2.getLocation().getMax() - referencePosition.getMin());
		int ref0 = Math.min(ref0Start, ref0End);

		int ref1Start, ref1End;
		ref1Start = Math.abs(o1.getLocation().getMin() - referencePosition.getMin());
		ref1End = Math.abs(o1.getLocation().getMax() - referencePosition.getMin());
		int ref1 = Math.min(ref1Start,ref1End);

		if (ref0 < ref1) {
			return 1;
		} else if (ref0 > ref1) {
			return -1;
		} 

		return 0;			
	}
}