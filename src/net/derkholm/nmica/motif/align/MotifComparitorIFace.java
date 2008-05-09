package net.derkholm.nmica.motif.align;

import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.motif.Motif;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;

public interface MotifComparitorIFace {
	public abstract Matrix2D bestHitsMatrix(Motif[] set0, Motif[] set1)
			throws Exception;

	public abstract Matrix2D bestHitsMatrix(Motif[] set) throws Exception;

	public abstract MotifPair[] allHits(Motif[] set0, Motif[] set1,
			double threshold) throws Exception;

	public abstract MotifPair[] bestHits(Motif[] set0, Motif[] set1)
			throws Exception;

	public abstract MotifPair[] bestReciprocalHits(Motif[] set0, Motif[] set1)
			throws Exception;

	public abstract double compareMotifs(WeightMatrix wm0, Distribution pad0,
			WeightMatrix wm1, Distribution pad1) throws Exception;

}