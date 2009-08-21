package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

import org.biojava.bio.BioError;
import org.bjv2.util.cli.Option;

abstract public class FilteringSAMProcessor extends SAMProcessor {
	protected File outFile = null;
	private String outString = "-";
	protected SAMFileWriter outWriter;
	protected boolean sorted = false;

	@Option(help="Output map file (file extension will decide if it's going to be written as SAM or BAM). If unspecified, or '-' specified, output made to stdout.", optional=true)
	public void setOut(String str) {
		this.outString = str;
	}
	
	@Option(help="Data is sorted (default = false)", optional=true)
	public void setSorted(boolean b) {
		this.sorted  = b;
	}
	
	public void initializeSAMWriter(boolean presorted) {
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		SAMFileHeader samHeader = new SAMFileHeader();
		
		if (inReader == null) {
			throw new BioError("Internal error: SAM reader needs to be initialised before the SAM writer");
		}
		
		if (this.outString.equals("-")) {
			this.outWriter = factory.makeSAMWriter(inReader.getFileHeader(), presorted, System.out);
		} else {
			this.outFile = new File(outString);
			this.outWriter = factory.makeSAMOrBAMWriter(inReader.getFileHeader(), presorted, this.outFile);
		}
	}
}