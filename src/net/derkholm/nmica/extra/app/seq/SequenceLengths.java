package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.NoSuchElementException;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Output the length of each input sequence", generateStub = true)
@NMExtraApp(launchName = "nmseqlen", vm = VirtualMachine.SERVER)
public class SequenceLengths {

	private RichSequenceIterator seqs;

	@Option(help="Input FASTA file")
	public void setSeqs(File f) throws FileNotFoundException {
		this.seqs = RichSequence.IOTools.readFastaDNA(
			new BufferedReader(new FileReader(f)), null);
	}
	
	public void main(String[] args) throws BioException {
		while (this.seqs.hasNext()) {
			Sequence s = this.seqs.nextSequence();
			System.out.printf("%s\t%s%n",s.getName(),s.length());
		}
	}
}
