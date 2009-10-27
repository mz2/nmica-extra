package net.derkholm.nmica.model.analysis;


public class MotifROCAUCSummary implements Comparable {
		private String motifName;
		private double auc;
		private double bootstrapFraction;

		public MotifROCAUCSummary(String motifName, double auc, double bootstrapFraction) {
			this.motifName = motifName;
			this.auc = auc;
			this.bootstrapFraction = bootstrapFraction;
		}

		public String getMotifName() {
			return motifName;
		}

		public void setMotifName(String motifName) {
			this.motifName = motifName;
		}

		public double getAuc() {
			return auc;
		}

		public void setAuc(double auc) {
			this.auc = auc;
		}

		public double getBootstrapFraction() {
			return bootstrapFraction;
		}

		public void setBootstrapFraction(double bootstrapFraction) {
			this.bootstrapFraction = bootstrapFraction;
		}

		public int compareTo(Object o) {
			if (!(o instanceof MotifROCAUCSummary)) {
				throw new IllegalArgumentException(
					"Trying to compare a MotifROCAUCSummary instance against a " + 
					o.getClass().getCanonicalName());
			}
			
			MotifROCAUCSummary sum = (MotifROCAUCSummary) o;
			return Double.compare(this.getAuc(), sum.getAuc());
		}
	}