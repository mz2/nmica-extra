package net.derkholm.nmica.extra.app;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Spike motifs to sequences with a specified rate", generateStub = true)
@NMExtraApp(launchName = "nmspikeseq", vm = VirtualMachine.SERVER)
public class Spike {
	protected File[] seqFiles;
	protected File[] motifFiles;
	protected String outFile;
	protected double rate;
	private Object type;

	@Option(help="Input sequence file(s)")
	public void setSeqs(File[] seqs) {
		seqFiles = seqs;
	}
	
	@Option(help="Input motif file(s)")
	public void setMotif(File[] motif) {
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
	
	@Option(help="Sequence type:dna(default)|protein")
	public void setType(String t) {
		type = t;
	}
	
	public void main(String[] args) throws Exception {
		FiniteAlphabet alp = ProteinTools.getAlphabet();
		
		Motif[] mot  = MotifIOTools.loadMotifSetXML(motifFiles);
		
		for (File seqFile : seqFiles) {
			BufferedReader seqReader = new BufferedReader(new FileReader(seqFile));
			SequenceIterator seqIt = SeqIOTools.readFastaProtein(seqReader);

			WeightMatrix wm1 =  mot[0].getWeightMatrix();

			//SequenceDB db = new HashSequenceDB();
			SequenceDB dbRandom = new HashSequenceDB();
			int count=0;
			int total=0;
			while (seqIt.hasNext())
				total++;
			
			seqReader = new BufferedReader(new FileReader(seqFile));
			if (type.equals("DNA"))
				seqIt = SeqIOTools.readFastaDNA(seqReader);
			else if (type.equals("protein"))
				seqIt = SeqIOTools.readFastaProtein(seqReader);
			else {
				System.err.println("Sequence types allowed: dna|protein");
				System.exit(1);
			}
			List<Sequence> idList = new ArrayList<Sequence>();
			
			while (seqIt.hasNext()) {
				Sequence seq = seqIt.nextSequence();
			        Sequence spikedSeq = seq;
				  count++;
				  double fraction = count * 1.0 / total; 
				  if (fraction<=rate) {
				      SymbolList smallMotifLikeSeq = generateSeqFromWM(wm1);
				      spikedSeq = insertSeq(smallMotifLikeSeq,seq,alp);
				      idList.add(spikedSeq);
				  } else break;
			}
			
			for (Sequence s : idList)
				dbRandom.addSequence(s);
			
			OutputStream os = new FileOutputStream(outFile);
			OutputStream outSeqsBuffer = new BufferedOutputStream(os);
			SeqIOTools.writeFasta(outSeqsBuffer,dbRandom);
			
			outSeqsBuffer.close();
		}		
	}

	private static Sequence insertSeq (SymbolList shortSeq, Sequence bigSeq, FiniteAlphabet alp) throws IllegalSymbolException{
		Random generator = new Random();
		
		int pos = generator.nextInt(bigSeq.length()-shortSeq.length());
		List<Symbol> l = new ArrayList<Symbol>();
		for (int i = 1; i <=pos;i++) {
			l.add(bigSeq.symbolAt(i));
		}
		l.addAll(shortSeq.toList());
		for (int i=pos+1;i<bigSeq.length();i++) {
			l.add(bigSeq.symbolAt(i));
		}
		return (new SimpleSequence(new SimpleSymbolList(alp,l),"",bigSeq.getName(),Annotation.EMPTY_ANNOTATION ));
	}
	
	private static SymbolList generateSeqFromWM (WeightMatrix wm) throws IllegalSymbolException{
		Random generator = new Random();
		Map<Symbol,Double> map = new HashMap<Symbol,Double>();
		Symbol[] spiky = new Symbol[wm.columns()];
		for (int pos = 0; pos<wm.columns();++pos) {
			Distribution d = wm.getColumn(pos);
			spiky[pos]=d.sampleSymbol();
			
		}
		SimpleSymbolList spikeWord = new SimpleSymbolList(spiky, spiky.length,ProteinTools.getAlphabet());
		
		return spikeWord;
	}
}

