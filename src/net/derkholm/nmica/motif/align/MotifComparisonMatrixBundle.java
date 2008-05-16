package net.derkholm.nmica.motif.align;

import net.derkholm.nmica.matrix.IntMatrix2D;
import net.derkholm.nmica.matrix.Matrix2D;

 public class MotifComparisonMatrixBundle {
    	public MotifComparisonMatrixBundle(
    			Matrix2D matrix, Matrix2D matrix2,
				IntMatrix2D offsetMatrix, IntMatrix2D offsetMatrix2) {
			this.senseScoreMatrix = matrix;
			this.antisenseScoreMatrix = matrix2;
			this.senseOffsetMatrix = offsetMatrix;
			this.antisenseOffsetMatrix = offsetMatrix2;
		}
    	
		protected Matrix2D senseScoreMatrix;
    	protected Matrix2D antisenseScoreMatrix;
    	protected IntMatrix2D senseOffsetMatrix;
    	protected IntMatrix2D antisenseOffsetMatrix;
    	
		public Matrix2D getSenseScoreMatrix() {
			return senseScoreMatrix;
		}
		public void setSenseScoreMatrix(Matrix2D senseScoreMatrix) {
			this.senseScoreMatrix = senseScoreMatrix;
		}
		public Matrix2D getAntisenseScoreMatrix() {
			return antisenseScoreMatrix;
		}
		public void setAntisenseScoreMatrix(Matrix2D antisenseScoreMatrix) {
			this.antisenseScoreMatrix = antisenseScoreMatrix;
		}
		public IntMatrix2D getSenseOffsetMatrix() {
			return senseOffsetMatrix;
		}
		public void setSenseOffsetMatrix(IntMatrix2D senseOffsetMatrix) {
			this.senseOffsetMatrix = senseOffsetMatrix;
		}
		public IntMatrix2D getAntisenseOffsetMatrix() {
			return antisenseOffsetMatrix;
		}
		public void setAntisenseOffsetMatrix(IntMatrix2D antisenseOffsetMatrix) {
			this.antisenseOffsetMatrix = antisenseOffsetMatrix;
		}
    }