package net.derkholm.nmica.extra.app.motif;

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

import org.biojava.bio.BioError;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Output XMS files in other formats (input given with unnamed arguments)", generateStub = true)
@NMExtraApp(launchName = "nmconvertxms", vm = VirtualMachine.SERVER)
public class MotifSetConverter {
	private OutputStream outputStream = System.out;
	
	private Format inputFormat = Format.XMS;
	private Format outputFormat = Format.TRANSFAC;
	
	private enum Format {
		XMS,
		TRANSFAC,
		MEME
	}
	
	@Option(help="Input format (default = xms)",optional=true)
	public void setInputFormat(Format format) {
		this.inputFormat = format;
	}

	@Option(help = "Output filename", optional = true)
	public void setOut(String str) throws FileNotFoundException {
		this.outputStream = 
			new BufferedOutputStream(new FileOutputStream(new File(str)));
	}
	
	@Option(help = "Output format (default = transfac)", 
			optional = true)
	public void setOutputFormat(Format format) {
		this.outputFormat = format;
	}

	public void main(String[] args) throws Exception {
		List<Motif> motifList = new ArrayList<Motif>();
		for (String filen : args) {
			FileReader fileReader = new FileReader(new File(filen));
			
			Motif[] motifs;
			if (inputFormat == Format.XMS) {
				motifs = MotifIOTools.loadMotifSetXML(fileReader);
			} else if (inputFormat == Format.MEME) {
				motifs = new Motif[]{MotifIOTools.loadMotifMEME(fileReader)};
			} else {
				motifs = null;
				throw new BioError("Unsupported format " + inputFormat);
			}
			
			if (fileReader != null)
				fileReader.close();
			for (Motif m : motifs)
				motifList.add(m);				
		}

		if (outputFormat == Format.TRANSFAC) {
			MotifIOTools.writeMotifSetTRANSFACFormat(
					outputStream, 
					motifList.toArray(new Motif[motifList.size()]));			
		} else {
			MotifIOTools.writeMotifSetXML(outputStream, 
					motifList.toArray(new Motif[motifList.size()]));
		}
	}
}
