package net.derkholm.nmica.extra.app;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;



@App(overview="A tool for regularising XMS motif set files", 
		generateStub=true)
@NMExtraApp(launchName="nmregularise", vm=VirtualMachine.SERVER)
public class MotifSetRegulariser {
	
	private String out;
	private double pseudoCount;
	private boolean ignoreErrors;
	
	
	@Option(help="Output filename")
	public void setOut(String str) {
		this.out = str;
	}
	
	@Option(help="Pseudocount (default=0.050)", optional=true)
	public void setPseudo(double d) {
		this.pseudoCount = d;
	}
	
	
	public void main(String[] args) throws Exception {
		List<Motif> allMotifs = new ArrayList<Motif>();

		for (String filen : args)
			for (Motif m : MotifIOTools.loadMotifSetXML(
							new BufferedReader(new FileReader(filen)))) 
				allMotifs.add(m);
		
		for (Motif m : allMotifs) {
			addPseudoCounts(m.getWeightMatrix(), pseudoCount);
			rescaleWeightMatrix(m.getWeightMatrix());
		}
		
		Motif[] motifs = allMotifs.toArray(new Motif[allMotifs.size()]);
		if (out != null)
			MotifIOTools.writeMotifSetXML(new BufferedOutputStream(
					new FileOutputStream(out)), motifs);
		else
			MotifIOTools.writeMotifSetXML(System.out, motifs);
	}
	
	public static void addPseudoCounts(WeightMatrix matrix, double pseudoCount) 
	throws IllegalSymbolException, ChangeVetoException {
		for (int i = 0; i < matrix.columns(); i++) {
			Distribution distrib = matrix.getColumn(i);
			for (Iterator j = ((FiniteAlphabet) distrib.getAlphabet()).iterator(); 
				j.hasNext(); ) {
				Symbol s = (Symbol)j.next();
				distrib.setWeight(s,distrib.getWeight(s) + pseudoCount);
				}
			}
			
		rescaleWeightMatrix(matrix);
		return;
	}

	public static void rescaleWeightMatrix(WeightMatrix matrix) throws IllegalSymbolException {
		for (int i = 0; i < matrix.columns(); i++) {
			Distribution distrib = matrix.getColumn(i);
			double sum = 0;
			for (Iterator j = ((FiniteAlphabet) distrib.getAlphabet()).iterator(); 
					j.hasNext(); )
				sum = sum + distrib.getWeight((Symbol)j.next());
			for (Iterator j = ((FiniteAlphabet) distrib.getAlphabet()).iterator(); 
					j.hasNext(); ) {
				Symbol s = (Symbol)j.next();
				if (sum > 0)
					distrib.setWeight(s,distrib.getWeight(s) / sum);
				else 
					throw new IllegalArgumentException(
						"Sum of the relative frequencies in the distribution <= 0");
			}
		}
		return;
	}

}
