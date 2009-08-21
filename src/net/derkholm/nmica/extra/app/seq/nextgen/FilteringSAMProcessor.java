package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;

import org.bjv2.util.cli.Option;

abstract public class FilteringSAMProcessor extends SAMProcessor {
	private File outFile = null;

	@Option(help="Output map file (file extension will decide if it's going to be written as SAM or BAM)")
	public void setOut(File f) {
		this.outFile  = f;
	}
}
