package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import net.derkholm.nmica.build.NMApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@App(overview="A tool for outputting subsets of motif sets based on name lists or a regular expression", generateStub=true)
@NMApp(launchName="nmgrep", vm=VirtualMachine.SERVER)
public class MotifGrep {

	private File motifs;
	private File list;
	private Pattern pattern;
	private String[] names;
	
	@Option(help="...",optional=true)
	public void setList(File list) {
		this.list = list;
	}
	
	@Option(help="Allowed motif names",optional=true)
	public void setNames(String[] names) {
		this.names = names;
	}


	@Option(help="Regular expression to match against the motif names",optional=true)
	public void setExp(String str) {
		this.pattern = Pattern.compile(str);
	}
	
	@Option(help="...")
	public void setMotifs(File motifs) {
		this.motifs = motifs;
	}

	/**
	 * @param args
	 */
	public void main(String[] args) 
		throws Exception
	{
		Motif[] motifs = MotifIOTools.loadMotifSetXML(new FileReader(this.motifs));
		List<Motif> om = new ArrayList<Motif>();
		
		if (names != null && list != null) {
			System.err.println("Supply only either -names or -list");
			System.exit(1);
		}
		
		if (names != null) {
			for (String str : names) {
				for (Motif m : motifs) {
					if (m.getName().equals(str)) {
						om.add(m);
						break;
					}	
				}
			}
		} else if (list != null) {
			BufferedReader br = new BufferedReader(new FileReader(list));
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				for (Motif m : motifs) {
					if (m.getName().equals(line)) {
						om.add(m);
						break;
					}
				}
			}
		} else if (pattern != null) {
			for (Motif m : motifs) {
				if (pattern.matcher(m.getName()).find()) {
					om.add(m);
				}
			}
		} else {
			System.err.println("Need to supply either -list, -names or -exp");
			System.exit(1);
		}
		MotifIOTools.writeMotifSetXML(System.out, om.toArray(new Motif[0]));
	}

}
