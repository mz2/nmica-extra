package net.derkholm.nmica.extra.app;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Spike motifs to sequences with a specified rate", generateStub = true)
@NMExtraApp(launchName = "nmspikeseq", vm = VirtualMachine.SERVER)
public class MotifSpike {
	protected InputStream[] seqFiles;
	protected InputStream[] motifFiles;
	protected String outFile;
	protected double rate = 1.0;
	protected String type = "DNA";

	@Option(help="Input sequence file(s)")
	public void setSeqs(InputStream[] seqs) {
		seqFiles = seqs;
	}
	
	@Option(help="Input motif file(s)")
	public void setMotif(InputStream[] motif) {
		motifFiles = motif;
	}
	
	@Option(help="Output filename",optional=true)
	public void setOut(String str) {
		outFile = str;
	}
	
	@Option(help="Spike rate")
	public void setRate(double d) {
		rate = d;
	}
	
	@Option(help="Sequence type:dna(default)|protein",optional=true)
	public void setType(String t) {
		type = t;
	}
	
	public void main(String[] args) throws Exception {
		FiniteAlphabet alp = null;
		Motif[] mot  = MotifIOTools.loadMotifSetXML(motifFiles);
		List<SymbolList> allSymLists = new ArrayList<SymbolList>();
		
		for (InputStream seqStream : seqFiles) {
			SequenceDB seqDB;
			if (type.equals("DNA")) {
				alp = DNATools.getDNA();
			}
			else if (type.equals("protein")) {
				alp = ProteinTools.getAlphabet();
			}
			else {
				System.err.println("Invalid sequence type. Types allowed:dna|protein");
				System.exit(1);
			}
			seqDB = SeqIOTools.readFasta(seqStream, alp);
			SequenceIterator seqIterator = seqDB.sequenceIterator();
			
			while (seqIterator.hasNext()) {
				allSymLists.add(new SimpleSymbolList(seqIterator.nextSequence()));
			}
		}
		
		SequenceDB spikedSeqDB = new HashSequenceDB();
		for (int i = 0; i < allSymLists.size(); i++) {
			SymbolList symList = allSymLists.get(i);
			spikedSeqDB.addSequence(
				new SimpleSequence(
					symList, null, "seq"+i, null));
		}
		
		for (InputStream inputStream : motifFiles) {
			WeightMatrix wm1 =  mot[0].getWeightMatrix();

			int seqCount = allSymLists.size();
			int spikeCount = (int) Math.round(rate * seqCount);
			
			Random random = new Random();
			while (spikeCount > 0) {
				int randIndex = random.nextInt(allSymLists.size());
			    SymbolList seq = allSymLists.get(randIndex);
				
			    SymbolList smallMotifLikeSeq = generateSeqFromWM(wm1);
			    insertSeqRandomlyToSeq(smallMotifLikeSeq, seq, alp);
			    spikeCount--;
			}
			
			OutputStream output;
			if (outFile == null) {
				output = System.out;
			} else {
				output = new BufferedOutputStream(new FileOutputStream(outFile));
			}
			SeqIOTools.writeFasta(output,spikedSeqDB);			
			
		}
	}

	private static void insertSeqRandomlyToSeq (SymbolList shortSeq, SymbolList seq, FiniteAlphabet alp) 
		throws 	IllegalSymbolException, 
				IndexOutOfBoundsException, 
				IllegalAlphabetException, 
				ChangeVetoException {
		Random generator = new Random();
		int pos = generator.nextInt(seq.length()-shortSeq.length()-1) + 1;
		seq.edit(new Edit(pos, shortSeq.length(), shortSeq));
	}
	
	private static SymbolList generateSeqFromWM (WeightMatrix wm) throws IllegalSymbolException {
		Random generator = new Random();
		Map<Symbol,Double> map = new HashMap<Symbol,Double>();
		Symbol[] spiky = new Symbol[wm.columns()];
		for (int pos = 0; pos<wm.columns();++pos) {
			Distribution d = wm.getColumn(pos);
			spiky[pos]=d.sampleSymbol();
			
		}
		SimpleSymbolList spikeWord = 
			new SimpleSymbolList(spiky, spiky.length,wm.getAlphabet());
		
		return spikeWord;
	}
}

