package net.derkholm.nmica.extra.app;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="A tool for separating XMS motif set files", 
		generateStub=true)
@NMExtraApp(launchName="nmseparate", vm=VirtualMachine.SERVER)
public class MotifSetSeparator {
	
	private String prefix;
	private boolean ignoreErrors;
	
	@Option(help="Output filename prefix",optional=true)
	public void setPrefix(String str) {
		this.prefix = str;
	}
	
	@Option(help="Ignore errors, e.g. empty/incorrect motif set files (default=false)",
			optional=true)
	public void setIgnoreErrors(boolean b) {
		this.ignoreErrors = b;
	}
	
	public void main(String[] args) throws Exception {
		List<Motif> motifList = new ArrayList<Motif>();
		for (String filen : args) {
			try {
				FileReader fileReader = new FileReader(new File(filen));
				Motif[] motifs = MotifIOTools.loadMotifSetXML(fileReader);
				if (fileReader != null) fileReader.close();
				for (Motif m : motifs) motifList.add(m);
			} catch (Exception e) {
				System.err.println(
					"ERROR! Could not parse motifs from " + filen);
				if (!ignoreErrors) {
					throw e;
				}
			}
		}
		
		for (Motif m : motifList) {
			String filen;
			if (prefix != null)
				filen = prefix + m.getName();
			else
				filen = m.getName();
			
			MotifIOTools.writeMotifSetXML(
				new BufferedOutputStream(
					new FileOutputStream(
						new File(filen))), motifList.toArray(
							new Motif[] {m}));
		}
	}
}
