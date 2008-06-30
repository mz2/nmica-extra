package net.derkholm.nmica.motif.align;

import net.derkholm.nmica.matrix.IntMatrix2D;
import net.derkholm.nmica.matrix.Matrix2D;
import net.derkholm.nmica.model.BitMatrix2D;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifComparisonMatrixBundle;

@Deprecated
public class SparseMotifComparisonMatrixBundle extends
		MotifComparisonMatrixBundle {
	protected BitMatrix2D allowedCoords;
	
	public SparseMotifComparisonMatrixBundle(Matrix2D matrix, Matrix2D matrix2,
			IntMatrix2D offsetMatrix, IntMatrix2D offsetMatrix2, 
			Motif[] rowMotifs, Motif[] colMotifs) {
		super(matrix, matrix2, offsetMatrix, offsetMatrix2, rowMotifs, colMotifs);
		
		allowedCoords = new BitMatrix2D(matrix.rows(), matrix.columns());
		for (int i = 0, rows = allowedCoords.rows(); i < rows; i++)
			for (int j = 0, cols = allowedCoords.columns(); j < cols; j++)
				allowedCoords.setValue(i, j, true);
	}
	
	public SparseMotifComparisonMatrixBundle(
			MotifComparisonMatrixBundle comparisonMatrixWithOffsets) {
		this(comparisonMatrixWithOffsets.getSenseScoreMatrix(),
			 comparisonMatrixWithOffsets.getAntisenseScoreMatrix(),
			 comparisonMatrixWithOffsets.getSenseOffsetMatrix(),
			 comparisonMatrixWithOffsets.getAntisenseOffsetMatrix(),
			 comparisonMatrixWithOffsets.getRowMotifs(),
			 comparisonMatrixWithOffsets.getColumnMotifs());
	}

	protected boolean isAllowed(int row, int column) {
		return allowedCoords.getValue(row, column);
	}
	
	protected void setAllowed(int row, int column, boolean v) {
		if (!v) {
			System.out.println("set row " + row + " and " + column + " not allowed.");
		}
		allowedCoords.setValue(row, column, v);
	}

}