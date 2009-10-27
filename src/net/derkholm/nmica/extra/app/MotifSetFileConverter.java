package net.derkholm.nmica.extra.app;

import java.io.File;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Convert motifset files in between formats", generateStub=true)
@NMExtraApp(launchName="nmconvert", vm=VirtualMachine.SERVER)
public class MotifSetFileConverter {

	String fromFormat = "transfac";
	String toFormat = "xms";
	String outFile = null;
	private File matrixF;
	private File factorF;
	private File classF;
	private File classificationF;
	private boolean outputUnclassified;
	private String nameFormat;
	
	private static final String MNAME_FORMAT_UNMODIFIED = "unmodified";
	private static final String MNAME_FORMAT_NAME_CLASS = "name_class";
	private static final String MNAME_FORMAT_NAME_CLASS_SPECIES = "name_class_species";
	
	@Option(help="The source format (default=transfac)",optional=true)
	public void setFrom(String str) {
		this.fromFormat = str;
	}
	
	@Option(help="The target format (default=xms)",optional=true)
	public void setTo(String str) {
		this.toFormat = str;
	}
	
	@Option(help="Input file(s) -- in case of TRANSFAC matrix parsing, this should be matrix.dat")
	public void setIn(File[] f) {
		
	}
	
	@Option(help="Output file")
	public void setOut(String fStr) {
		this.outFile = fStr;
	}
	
	@Option(help="Path to TRANSFAC factor.dat",optional=true)
	public void setFactors(File f) {
		this.factorF = f;
	}
	
	@Option(help="Path to TRANSFAC class.dat",optional=true)
	public void setClasses(File f) {
		this.classF = f;
	}
	
	@Option(help="Path to TRANSFAC classification.dat",optional=true)
	public void setClassifications(File f) {
		this.classificationF = f;
	}
	
	@Option(help="Output also motifs with no TRANSFAC classification",optional=true)
	public void setOutputUnclassified(boolean b) {
		this.outputUnclassified = b;
	}
	
	@Option(help="Motif name format:unmodified|name_class|name_class_species (default=unmodified)")
	public void setNameFormat(String str) {
		this.nameFormat = str;
	}
	
	public void main(String[] argv) {
		
	}
}
