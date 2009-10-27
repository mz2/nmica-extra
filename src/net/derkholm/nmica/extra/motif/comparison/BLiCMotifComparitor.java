package net.derkholm.nmica.extra.motif.comparison;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;

import net.derkholm.nmica.motif.MotifComparitor;
import net.derkholm.nmica.motif.MotifComparitorIFace;

public class BLiCMotifComparitor extends MotifComparitor {
	private static BLiCMotifComparitor instance;
	
	@Override
	public double compareMotifs(WeightMatrix wm0, Distribution pad0,
			WeightMatrix wm1, Distribution pad1) {
		return 0;
		
	}

	@Override
	public ScoreOffsetPair compareMotifsWithOffset(WeightMatrix wm0,
			Distribution pad0, WeightMatrix wm1, Distribution pad1) {
		return null;
	}

	public static MotifComparitorIFace getMotifComparitor() {
		if (instance == null)
			instance = new BLiCMotifComparitor();
		
		return instance;
	}

}
