/**
 * 
 */
package net.derkholm.nmica.extra.peak;

import java.util.Comparator;

import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakEntry;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakEntryComparator;
import net.derkholm.nmica.extra.app.seq.nextgen.RetrievePeakSequencesFromEnsembl.PeakFormat;

public class PeakEntryAscComparitor extends PeakEntryComparator implements Comparator<PeakEntry> {
	PeakFormat format;
	PeakEntryDescComparitor descComparitor;
	private boolean seqName;
	
	public PeakEntryAscComparitor(boolean includeSeqName, PeakFormat format) {
		this.seqName = includeSeqName;
		this.format = format;
		this.descComparitor = new PeakEntryDescComparitor(includeSeqName,format);
	}

	public int compare(PeakEntry o1, PeakEntry o2) {
		return -descComparitor.compare(o1, o2);
	}
}