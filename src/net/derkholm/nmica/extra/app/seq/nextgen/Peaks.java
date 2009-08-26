package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@NMExtraApp(launchName = "ngpeaks", vm = VirtualMachine.SERVER)
@App(overview = "Call peak regions from sequencing data", generateStub = true)
public class Peaks extends SAMProcessor {

	private File depthFile;
	private File controlDepthFile;
	private double pValue;
	private double foldChange;

	@Option(help="Sequencing depth of the ChIP reads (ngdepth file)")
	public void setDepth(File f) {
		this.depthFile = f;
	}
	
	@Option(help="Sequencing depth of the control reads (ngdepth file)")
	public void setControlDepth(File f) {
		this.controlDepthFile = f;
	}
	
	@Override
	@Option(help="Target window size for ChIP reads and the control reads")
	public void setWindowSize(int i) {
		super.setWindowSize(i);
	}
	
	@Option(help="P-value cutoff for reporting peaks")
	public void setPValue(double d) {
		this.pValue = d;
	}
	
	@Option(help="Fold change cutoff")
	public void setFoldChange(double d) {
		this.foldChange = d;
	}
	
	public void main(String[] args) {
		
	}
}