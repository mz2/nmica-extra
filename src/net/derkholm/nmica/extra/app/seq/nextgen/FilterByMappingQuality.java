package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Filter mapped reads by mapping quality", generateStub = true)
@NMExtraApp(launchName = "ngfiltermap", vm = VirtualMachine.SERVER)
public class FilterByMappingQuality {

	private SAMFileReader sam;
	private int mapQual = 10;
	private File outFile;


	@Option(help="Sequencing reads mapped to reference in SAM/BAM format. " +
			"If '-' specified, will read from stdin.")
	public void setMap(String str) {
		if (str.equals("-")) {
			this.sam = new SAMFileReader(System.in);
		} else {
			this.sam = new SAMFileReader(new File(str));
		}
		sam.setValidationStringency(ValidationStringency.SILENT);
	}
	
	@Option(help="Mapping quality threshold (default = 10)")
	public void setMappingQualityAbove(int i) {
		this.mapQual = i;
	}
	

	
	public void main(String[] args) throws FileNotFoundException {
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		SAMFileWriter writer = factory.makeSAMOrBAMWriter(this.sam.getFileHeader(), false, this.outFile);
		
		OutputStream stream = null;
		if (outFile == null) {
			stream = System.out;
		} else {
			stream = new FileOutputStream(outFile);
		}
		for (final SAMRecord rec : this.sam) {
			if (rec.getMappingQuality() < this.mapQual) {
				writer.addAlignment(rec);
			}
		}
		writer.close();
	}
}
