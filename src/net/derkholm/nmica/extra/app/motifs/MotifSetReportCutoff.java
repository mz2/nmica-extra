package net.derkholm.nmica.extra.app.motifs;

import java.io.FileNotFoundException;
import java.io.FileReader;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.bjv2.util.cli.App;

@App(overview="Report the score threshold annotated for each of the motifs in the input file(s)", generateStub=true)
@NMExtraApp(launchName="nmreportcutoff", vm=VirtualMachine.SERVER)
public class MotifSetReportCutoff {

	public void main(String[] args) throws FileNotFoundException, Exception {
		for (String filename : args) {
			Motif[] motifs = MotifIOTools.loadMotifSetXML(new FileReader(filename));
			
			for (Motif m : motifs) {
				System.out.printf("%s\t%s%n",m.getName(),m.getThreshold());
			}
		}
	}
}