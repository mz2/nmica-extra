package net.derkholm.nmica.extra.app;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStream;
import java.util.ArrayList;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.motif.NMWeightMatrix;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "A tool for concatenating columns from two or more motifs together", generateStub = true)
@NMExtraApp(launchName = "nmconcatenate", vm = VirtualMachine.SERVER)

public class MotifConcatenator {
	private File outputFile;
	private String name;

	@Option(help="Output file (written to stdout if not given)", optional=true)
	public void setOut(File of) {
		this.outputFile = of;
	}
	
	@Option(help="Output motif name (concatenation of motif names if not given)", optional=true)
	public void setName(String str) {
		this.name = str;
	}
	
	void main(String[] args) throws FileNotFoundException, Exception {
		ArrayList<Motif> motifs = new ArrayList<Motif>();
		int totalCount = 0;
		
		for (String filen : args) {
			File motiff = new File(filen);
			Motif[] ms = MotifIOTools.loadMotifSetXML(new FileReader(motiff));
			for (Motif m : ms) {
				motifs.add(m);
				totalCount = totalCount + m.getWeightMatrix().columns();
			}
		}
		
		Distribution[] dists = new Distribution[totalCount];
		{ 
			int i = 0;
			for (Motif m : motifs) {
				WeightMatrix wm = m.getWeightMatrix();
				
				for (int c = 0,cnt = wm.columns(); c < cnt; c++) {
					dists[i] = wm.getColumn(c);
					i++;
				}
			}
			assert i == (totalCount-1);
		}
		
		WeightMatrix wm = new NMWeightMatrix(dists, totalCount, 0);
		
		if (name == null) {
			StringBuffer nameBuf = new StringBuffer();
			
			for (int m=0,cnt=motifs.size() - 1; m < cnt; m++) {
				nameBuf.append(motifs.get(m).getName() + "_");
			}
			
			nameBuf.append(motifs.get(motifs.size()-1).getName());
			name = nameBuf.toString();
		}
		
		Motif m = new Motif();
		m.setName(name);
		m.setWeightMatrix(wm);
		
		OutputStream ostr;
		if (outputFile != null) {
			ostr = new FileOutputStream(outputFile);
		} else {
			ostr = System.out;
		}
		
		MotifIOTools.writeMotifSetXML(ostr, new Motif[] {m});
	}
}
