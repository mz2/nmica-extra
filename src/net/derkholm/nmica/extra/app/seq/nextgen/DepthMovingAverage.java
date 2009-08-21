package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.sf.samtools.SAMRecord;

import org.biojava.bio.BioException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

@NMExtraApp(launchName = "ngdepth", vm = VirtualMachine.SERVER)
@App(overview = "Output sequencing depth in a moving average", generateStub = true)
public class DepthMovingAverage extends SAMProcessor {
	
	public enum Format {
		SQLITE,
		TSV
	}

	
	
	private Format format = Format.TSV;

	private void setFormat(Format format) {
		this.format = format;
	}
	
	@Option(help="Query overlapping / contained reads (default=overlapping)", optional=true, userLevel=UserLevel.DEBUG)
	public void setQuery(QueryType queryType) {
		 super.setQueryType(queryType);
	}
	
	@Option(help="Iterate in moving window / overlapping reads", optional=true, userLevel=UserLevel.DEBUG)
	public void setIterate(IterationType iterType) {
		if (iterType == IterationType.ONE_BY_ONE) {
			System.err.println("Iteration type cannot be one_by_one");
			System.exit(3);
		}
		super.setIterationType(iterType);
	}
	
	public void main(String[] args) throws BioException {
		setIterationType(IterationType.MOVING_WINDOW);
		
		initializeSAMReader();
		process();
	}
	
	@Override
	public void process(final List<SAMRecord> recs, String refName, int begin, int end, int seqLength) {
		double avg = 0.0;
		
		for (SAMRecord rec : recs) {
			
			if (format == Format.TSV) {
				
			}
			
		}
	}
}