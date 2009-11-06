//
//  MotifPair.java
//  MotifExplorer
//
//  Created by Thomas Down on 08/03/2005.
//  Copyright 2005 __MyCompanyName__. All rights reserved.
//

package net.derkholm.nmica.motif.align;

import net.derkholm.nmica.motif.Motif;

public class MotifPair implements Comparable {
	    protected final Motif m1;
	    protected final Motif m2;
	    protected final boolean flip;
	    protected final double score;
	    protected final double pValue;
	    
	    public MotifPair(Motif m1, Motif m2, double score, double pValue, boolean flip) {
	        this.m1 = m1;
	        this.m2 = m2;
	        this.score = score;
	        this.flip = flip;
	        this.pValue = pValue;
	    }

		public int compareTo(Object o) {
            MotifPair po = (MotifPair) o;
            double dif = score - po.score;
            if (dif < 0) {
                return 1;
            } else if (dif > 0) {
                return -1;
            } else {
                return 0;
            }
        }
        
        public boolean equals(Object o) {
            MotifPair po = (MotifPair) o;
            return m1 == po.m1 && m2 == po.m2;
        }
        
        public int hashCode() {
            return m1.hashCode() * 37 + m2.hashCode();
        }


		public Motif getM1() {
			return m1;
		}


		public Motif getM2() {
			return m2;
		}


		public double getScore() {
			return score;
		}
	
		public boolean isFlipped() {
			return flip;
		}
}