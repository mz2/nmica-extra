package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import org.biojava.bio.BioException;
import org.bjv2.util.cli.App;

@NMExtraApp(launchName = "ngdepth", vm = VirtualMachine.SERVER)
@App(overview = "Output sequencing depth inside a window.", generateStub = true)
public class DepthMovingAverage extends SAMProcessor {
	
	public enum Format {
		SQLITE,
		TSV
	}
	
	private Format format = Format.TSV;
	private int windowIndex;
	private Map<String, Integer> readCounts;

	private void setFormat(Format format) {
		this.format = format;
	}
	
	public void main(String[] args) throws BioException {
		setIterationType(IterationType.MOVING_WINDOW);
		setQueryType(QueryType.OVERLAP);
		initializeSAMReader();
		
		this.windowIndex = 0;
		process();
	}
	
	@Override
	public void process(final List<SAMRecord> recs, String refName, int begin, int end, int seqLength) {
		double avg = 0.0;
		int recCount = recs.size();
		if (recCount > 0) {
			System.out.printf("%s\t%d\t%d\t%d\t%d%n", refName, this.windowIndex, begin, end, recs.size());
		}
		
		this.windowIndex++;
	}
}