/**
 * 
 */
package net.derkholm.nmica.extra.peak;

import java.util.Comparator;

import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakEntry;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakEntryComparator;

public class PeakEntryRandomComparitor extends PeakEntryComparator implements Comparator<PeakEntry> {
	
	public PeakEntryRandomComparitor() {
		
	}
	
	public int compare(PeakEntry o1, PeakEntry o2) {
		return 0;
	}
}