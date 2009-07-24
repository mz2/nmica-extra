package net.derkholm.nmica.extra.app.motifs;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

@App(overview = "Output XMS files in other formats", generateStub = true)
@NMExtraApp(launchName = "nmconvertxms", vm = VirtualMachine.SERVER)
public class MotifSetConverter {
	private OutputStream outputStream = System.out;
	private OutputFormat outputFormat = OutputFormat.TRANSFAC;
	
	private enum OutputFormat {
		TRANSFAC
	}

	@Option(help = "Output filename", optional = true)
	public void setOut(String str) throws FileNotFoundException {
		this.outputStream = 
			new BufferedOutputStream(new FileOutputStream(new File(str)));
	}
	
	@Option(help = "Output format (default = TRANSFAC)", 
			userLevel = UserLevel.EXPERT, 
			optional = true)
	public void setFormat(OutputFormat format) {
		this.outputFormat = format;
	}

	public void main(String[] args) throws Exception {
		List<Motif> motifList = new ArrayList<Motif>();
		for (String filen : args) {
			FileReader fileReader = new FileReader(new File(filen));
			Motif[] motifs = MotifIOTools.loadMotifSetXML(fileReader);
			if (fileReader != null)
				fileReader.close();
			for (Motif m : motifs)
				motifList.add(m);
		}

		if (outputFormat == OutputFormat.TRANSFAC) {
			MotifIOTools.writeMotifSetTRANSFACFormat(
					outputStream, 
					motifList.toArray(new Motif[motifList.size()]));			
		}
	}
}
