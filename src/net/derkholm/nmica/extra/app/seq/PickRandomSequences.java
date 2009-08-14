package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Write out random sequences or sections of sequences", generateStub = true)
@NMExtraApp(launchName = "nmrandomseq", vm = VirtualMachine.SERVER)
public class PickRandomSequences {
	private int count = 1;
	private File seqFile;
	private int length = -1;
	private boolean sampleWithReplacement;
	private boolean sampleFragsWithReplacement;
	private static Random random = new Random();
	
	@Option(help="Number of sequences to sample (default=1)", optional=true)
	public void setCount(int count) {
		this.count = count;
	}
	
	@Option(help="Input sequences")
	public void setSequences(File seqFile) {
		this.seqFile = seqFile;
	}
	
	@Option(help="Fragment length to sample (default=-1, interpreted as the length of the sequence sampled)", optional=true)
	public void setLength(int length) {
		this.length = length;
	}
	
	@Option(help="Sample with replacement from sequences, i.e. allow same sequence in output (default=false)", optional=true)
	public void setSampleSeqsWithReplacement(boolean b) {
		this.sampleWithReplacement = b;
	}
	
	@Option(help="Sample with replacement from sequence fragments, i.e. allow the same position in input sequences more than once")
	public void setSampleFragsWithReplacement(boolean b) {
		this.sampleFragsWithReplacement = b;
	}
	
	public void main(String[] args) throws BioException, IOException {
		RichSequenceIterator seqIterator = RichSequence.IOTools.readFastaDNA(new BufferedReader(new FileReader(seqFile)), null);
		
		List<Sequence> seqs = new ArrayList<Sequence>();
		while(seqIterator.hasNext()) {
			seqs.add(seqIterator.nextSequence());
		}
		
		List<Sequence> chosenSeqs = new ArrayList<Sequence>();
		if (!sampleWithReplacement) {
			int c = count;			
			while (c > 0) {
				int randSeqIndex = random.nextInt(seqs.size());
				chosenSeqs.add(seqs.remove(randSeqIndex));
				
				c--;
			}
		} 
		else if (sampleWithReplacement || (length > 0)) {
			/* if you want to sample from sequences with replacement
			 * or if the wanted length is specified 
			 */
			
			if (length > 0) {
				
			}
			
		}
		
		for (Sequence seq : chosenSeqs) {
			RichSequence.IOTools.writeFasta(System.out, seq, null);
		}
		seqs = null;
	}
}