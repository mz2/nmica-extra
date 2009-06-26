package net.derkholm.nmica.extra.app;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.extra.motif.comparison.BLiCMotifComparitor;
import net.derkholm.nmica.extra.motif.comparison.KullbackLeiblerDifferenceMotifComparitor;
import net.derkholm.nmica.model.metamotif.Dirichlet;
import net.derkholm.nmica.model.metamotif.MetaMotif;
import net.derkholm.nmica.model.metamotif.MetaMotifIOTools;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifComparitorIFace;
import net.derkholm.nmica.motif.MotifIOTools;
import net.derkholm.nmica.motif.SquaredDifferenceMotifComparitor;
import net.derkholm.nmica.motif.align.InvalidMetaMotifException;
import net.derkholm.nmica.motif.align.MotifAlignment;

import org.biojava.bio.symbol.IllegalSymbolException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="A tool for aligning a set of motifs and output the alignment", 
	generateStub=true)
@NMExtraApp(launchName="nmalign", vm=VirtualMachine.SERVER)
public class MotifAligner {
	
	private Motif[] motifs;
	private String outputType = "align_cons";
	
	private String outFile;
	private PrintStream outputStream = System.out;
	
	private MotifComparitorIFace mc = SquaredDifferenceMotifComparitor.getMotifComparitor();
	private int minColPerPos = 2;
	private boolean outputSingleMotif;
	private double singleMotifPseudoCount;
	private double singleMotifPrecision = 10.0;
	private boolean addName;
	private double minColWeight;
	private double maximumPrecision;
	private String prefix;
	
	/*
	@Option(help="Output file",optional=true)
	public void setOut(String outFile) {
		this.outFile = outFile;
	}*/
	
	@Option(help="Add the name of the motifset file to the motif sets to make them unique (default=false)",
			optional=true)
	public void setAddName(boolean b) {
		this.addName = b;
	}
	
	@Option(help="Output file", optional=true)
	public void setOut(String str) throws FileNotFoundException {
		this.outFile = str;
		this.outputStream = new PrintStream(new FileOutputStream(new File(str)));
	}
	
	@Option(help="Add a freeform prefix to motif names",
			optional=true)
	public void setNamePrefix(String str) {
		this.prefix = str;
	}
	
	@Option(help="Minimum number of columns per position " +
			"to allow it to make it to output (default=2)",optional=true)
	public void setMinCols(int i) {
		this.minColPerPos = i;
	}
	
	@Option(help="If only one motif is given as an input, output it as is " +
			"(metamotif output type gets special treatment if this is specified, " +
			"see -singleMotifPseudoCount and -singleMotifPrecision)",optional=true)
	public void setOutputSingleMotif(boolean b ) {
		this.outputSingleMotif = b;
	}
	
	@Option(help="In the metamotif output mode weight of each column will be capped to this value",
			optional=true)
	public void setMinColWeight(double d) {
		this.minColWeight = d;
	}
	
	@Option(help="In the metamotif output mode precision of each column will be capped to this value",
			optional=true)
	public void setMaxPrecision(double d) {
		this.maximumPrecision = d;
	}
	
	@Option(help="If the output contains only a single motif and the metamotif output type is specified, " +
			"-singleMotifPseudoCount specifies the pseudo count applied to each output Dirichlet column.",optional=true)
	public void setSingleMotifPseudoCount(double d) {
		this.singleMotifPseudoCount = d;
	}
	
	@Option(help="If the output contains only a single motif and the metamotif output type is specified, " +
			"-singleMotifPrecision specifies the precision applied to each output Dirichlet column.",optional=true)
	public void setSingleMotifPrecision(double d) {
		this.singleMotifPrecision  = d;
	}
	
	@Option(help="Distance metric:sqdiff2.5(default)|sqdiff|kd|blic",optional=true)
	public void setDist(String dist) {
		if (!dist.equals("sqdiff2.5") &!
				dist.equals("sqdiff") &!
				dist.equals("kd") &!
				dist.equals("blic"))
			return;
		
		if (dist.equals("sqdiff2.5"))
			mc = SquaredDifferenceMotifComparitor.getMotifComparitor();
		else if (dist.equals("sqdiff"))
			mc = SquaredDifferenceMotifComparitor.getCarthesianMotifComparitor();
		else if (dist.equals("kd"))
			mc = KullbackLeiblerDifferenceMotifComparitor.getMotifComparitor();
		else if (dist.equals("blic"))
			mc = BLiCMotifComparitor.getMotifComparitor();
		else {
			System.err.println("Incorrect distance metric type : " + dist);
			System.exit(1);
		}
	}
	
	@Option(help="Output type (default=align_cons): avg|metamotif|cons|align_cons",optional=true)
	public void setOutputType(String type) {
		if (!type.equals("avg") &! 
				type.equals("all") &!
				type.equals("metamotif") &! 
				type.equals("cons") &!
				type.equals("align_cons")) {
			System.err.println("Allowed output types: all|avg|metamotif|cons|align_cons (default=align_cons)");
			System.exit(1);
		}
		this.outputType = type;
	}
	
	public void main(String[] args) throws FileNotFoundException, Exception {
		List<Motif> motifList = new ArrayList<Motif>();
		for (String fStr : args) {
			BufferedInputStream f = new BufferedInputStream(new FileInputStream(fStr));
			Motif[] ms = MotifIOTools.loadMotifSetXML(f);
			
			for (Motif m : ms) {
				if (this.addName)
					m.setName(m.getName() + "_" + fStr);
				else if (this.prefix != null)
					m.setName(prefix + m.getName());
			}

			motifList.addAll(Arrays.asList(ms));
		}
		this.motifs = motifList.toArray(new Motif[motifList.size()]);
		
		if (motifs.length == 1) {
			System.err.println(
				"WARNING! The supplied motif set file contains only one motif.\n");
			
			//if (outputSingleMotif &! outputType.equals("metamotif")) {
			//	System.out.println();
			//}
			
			//System.exit(1);
		}
		
		//TODO: Fix the alignment so you don't have to do this three times...
		MotifAlignment alignment 
			= new MotifAlignment(motifs, mc);
		alignment = new MotifAlignment(alignment.motifs(), mc);
		alignment = new MotifAlignment(alignment.motifs(), mc);
		alignment = alignment.alignmentWithZeroOffset();
		//alignment = alignment.trimToColumnsPerPosition(minColPerPos);
		if (minColPerPos > 1 && motifs.length > 1) {
			alignment = alignment.trimToColumnsPerPosition(minColPerPos);
		}
		
		//alignment.setName(this.motifs.getName().replace("\\.xms", "")+"_aligned");
				
		if (!outputType.equals("align_cons")) {
			System.err.println(alignment.alignmentConsensusString());
		}
		
		if (outputType.equals("avg"))
			MotifIOTools.writeMotifSetXML(
					outputStream,
					new Motif[] {alignment.averageMotif()});
		if (outputType.equals("all"))
			MotifIOTools.writeMotifSetXML(
					outputStream, 
					alignment.motifs());
		else if (outputType.equals("metamotif")) {
			//System.err.println(alignment.alignmentConsensusString());
			try {
				MetaMotif mm = alignment.metamotif(true);
				if (this.minColWeight > 0) MotifSetRegulariser.addPseudoCounts(mm, minColWeight);
				if (this.maximumPrecision > 0) capPrecision(mm, this.maximumPrecision);
				
				MetaMotifIOTools.
					writeMetaMotifSetToMotifSetWithAnnotations(
						outputStream, 
						new MetaMotif[] {mm});
			} catch (InvalidMetaMotifException e) {
				System.err.printf("Will output the input motif as a metamotif " +
				"with the per-column precision of %.3f.%n, (pseudocount=%.3f)",singleMotifPrecision, singleMotifPseudoCount);
				MetaMotifIOTools.
					writeMetaMotifSetToMotifSetWithAnnotations(
						outputStream, 
						new MetaMotif[] {
								new MetaMotif(motifs[0],
														this.singleMotifPrecision, 
														this.singleMotifPseudoCount)});				
			}
			
		} else if (outputType.equals("cons"))
			outputStream.println(alignment.consensusString());
		else if (outputType.equals("align_cons"))
			outputStream.println(alignment.alignmentConsensusString());

		System.err.println("Done.");
		outputStream.close();
		System.exit(0);
	}
	
	public static void capPrecision(MetaMotif mm, double maxPrecision) 
		throws IllegalSymbolException {
		System.err.printf("Setting precision to %.2f%n", maxPrecision);
		for (int i = 0; i < mm.columns(); i++) {
			Dirichlet distrib = mm.getColumn(i);
			if (distrib.alphaSum() > maxPrecision)
				distrib.setAlphaSum(maxPrecision);
		}
		//MotifSetRegulariser.rescaleWeightMatrix(mm);
		return;
	}
}