package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import net.derkholm.nmica.apps.MotifFinder;
import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.genematrix.GeneByGeneMatrix;
import net.derkholm.nmica.model.genematrix.GeneByGeneMatrixExpansionBundle;

import org.bjv2.util.cli.App;

@App(overview="Expand a gene-by-gene matrix", generateStub=true)
@NMExtraApp(launchName = "nmexpandgenematrix", vm = VirtualMachine.SERVER)
public class ExpandGeneByGeneMatrix {
	
	protected File matrixFile;
	protected File seqCountsFile;
	protected int[] seqCounts;

	public void setMatrix(File f) {
		this.matrixFile = f;
	}
	
	public void setSeqCountsFile(File f) {
		this.seqCountsFile = f;
	}
	
	public void main(String[] args) throws IOException {
		
		GeneByGeneMatrix m = MotifFinder.loadGeneByGeneMatrix(this.matrixFile, null);
		this.seqCounts = ExpandGeneByGeneMatrix.loadSeqCounts(this.seqCountsFile);
			
		GeneByGeneMatrixExpansionBundle expandedMBundle = m.expand(this.seqCounts);
		System.err.println(expandedMBundle.expandedMatrix);
	}
	
	public static int[] loadSeqCounts(File seqCountsFile) throws IOException {
		List<Integer> seqCountsList = new ArrayList<Integer>();
		BufferedReader r = new BufferedReader(new FileReader(seqCountsFile));
		
		String line = r.readLine();
		
		StringTokenizer stringTok = new StringTokenizer(line,"\t");
		while (stringTok.hasMoreTokens()) {
			seqCountsList.add(Integer.parseInt(stringTok.nextToken()));
		}
		
		int[] seqCounts = new int[seqCountsList.size()];
		for (int i = 0; i < seqCounts.length; i++) {seqCounts[i] = seqCountsList.get(i);}
		
		return seqCounts;
	}
	
}
