package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

import org.bjv2.util.cli.Option;

abstract public class FilteringSAMProcessor extends SAMProcessor {
	protected File outFile = null;
	private String outString = "-";
	protected SAMFileWriter outWriter;

	@Option(help="Output map file (file extension will decide if it's going to be written as SAM or BAM). If unspecified, or '-' specified, output made to stdout.", optional=true)
	public void setOut(String str) {
		this.outString = str;
	}
	
	public void initializeSAMWriter(boolean presorted) {
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		if (this.outString.equals("-")) {
			this.outWriter = factory.makeSAMWriter(null, presorted, this.outFile);
		} else {
			this.outFile = new File(outString);
			this.outWriter = factory.makeSAMOrBAMWriter(null, presorted, this.outFile);
		}
	}
}