package net.derkholm.nmica.extra.motif.comparison;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;

import net.derkholm.nmica.motif.MotifComparitor;

public class KullbackLeiblerDifferenceMotifComparitor extends MotifComparitor {

	public static MotifComparitor getMotifComparitor() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double compareMotifs(WeightMatrix wm0, Distribution pad0,
			WeightMatrix wm1, Distribution pad1) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public ScoreOffsetPair compareMotifsWithOffset(WeightMatrix wm0,
			Distribution pad0, WeightMatrix wm1, Distribution pad1) {
		// TODO Auto-generated method stub
		return null;
	}

}
