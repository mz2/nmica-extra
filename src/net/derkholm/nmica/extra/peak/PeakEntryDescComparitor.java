/**
 * 
 */
package net.derkholm.nmica.extra.peak;

import java.util.Comparator;

import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakEntry;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakEntryComparator;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakFormat;

public class PeakEntryDescComparitor extends PeakEntryComparator implements Comparator<PeakEntry> {
	PeakFormat format;
	boolean seqName;
	
	public PeakEntryDescComparitor(boolean includeSeqName, PeakFormat format) {
		this.seqName = includeSeqName;
		this.format = format;
	}
	
	public int compare(PeakEntry o1, PeakEntry o2) {
		if (seqName){
			int seqNameCompare = o1.seqName.compareTo(o2.seqName);
			if (seqNameCompare != 0) return seqNameCompare;
		}
		
		if (this.format == PeakFormat.SWEMBL || this.format == PeakFormat.FINDPEAKS) {
			int valCompare = Double.compare(o1.score, o2.score); /* Both o1 and o2 should be NaN for Swembl so this should return 0 */
			if (valCompare != 0) return valCompare;

			int fdrCompare = Double.compare(o1.fdr, o2.fdr);
			if (fdrCompare != 0) return fdrCompare;
			
			int tagCountCompare = Double.compare(o1.tagCount, o2.tagCount);
			if (tagCountCompare != 0) return tagCountCompare;

		} else if (this.format == PeakFormat.MACS) {
			int fdrCompare = Double.compare(o1.fdr, o2.fdr);
			if (fdrCompare != 0) return fdrCompare;

			int scoreCompare = -Double.compare(o1.score, o2.score);
			if (scoreCompare != 0) return scoreCompare;
			
			int tagCountCompare = Double.compare(o1.tagCount, o2.tagCount);
			if (tagCountCompare != 0) return tagCountCompare;
			
		} 
		
		return o1.id.compareTo(o2.id);
	}
}