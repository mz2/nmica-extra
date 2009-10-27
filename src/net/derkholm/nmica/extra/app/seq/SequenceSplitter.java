package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

/**
 * @author thomas,matias
 */
@App(overview = 
	"A tool for splitting sequences to " +
	"specified minimum and maximum lengths", generateStub = true)
@NMExtraApp(launchName = "nmsplitseq", vm = VirtualMachine.SERVER)
public class SequenceSplitter {
    private int maxLength = 200;
    private int minLength = 50;
    private boolean addSuffix = true;
    
    @Option(help="Maximum length (default=200)", optional=true)
    public void setMaxLength(int i) {
        this.maxLength = i;
    }
    
    @Option(help="Minimum length (default=50). " +
    		"Sequences that are shorter than this after splitting " +
    		"will not be output.", optional=true)
    public void setMinLength(int i) {
    	this.minLength = i;
    }
    
    @Option(help="Add a suffix to the (default=true)", optional=true)
    public void setAddSuffix(boolean b) {
    	this.addSuffix = b;
    }
    
    public void main(String[] args)
    		throws Exception {
    	if (args.length > 0) {
	        for (int i = 0; i < args.length; ++i) {
	            File seqFile = new File(args[i]);
	            doSplit(RichSequence.IOTools.readFastaDNA(new BufferedReader(new FileReader(seqFile)), null));
	        }
    	} else {
    		doSplit(RichSequence.IOTools.readFastaDNA(new BufferedReader(new InputStreamReader(System.in)), null));
    	}
    }
    
    private void doSplit(SequenceIterator si)
    	throws Exception
    {
         while (si.hasNext()) {
             Sequence seq = si.nextSequence();
             if (seq.length() > maxLength) {
                 int numChunks = (int) Math.ceil((1.0 * seq.length()) / maxLength);
                 int chunkSize = (int) Math.ceil((1.0 * seq.length()) / numChunks);
                 int pos = 1;
                 while (pos < seq.length()) {
                     int maxPos = Math.min(seq.length(), pos + chunkSize - 1);
                     if ((maxPos - pos) >= minLength) {
                    	 RichSequence.IOTools.writeFasta(System.out, new SimpleSequence(
                                 seq.subList(pos, maxPos),
                                 null,
                                 seq.getName() + "_" + pos + "-" + maxPos,
                                 Annotation.EMPTY_ANNOTATION), null);
                     }
                     pos = maxPos + 1;
                 }
             } else if (seq.length() > minLength) {
            	 RichSequence.IOTools.writeFasta(System.out, seq, null);
             }
         }
    }
    
    public static List<Sequence> splitSequences(SequenceIterator si, int minLength, int maxLength) throws Exception {
    	List<Sequence> seqs = new ArrayList<Sequence>();
    	
	     while (si.hasNext()) {
	         Sequence seq = si.nextSequence();
	         if (seq.length() > maxLength) {
	             int numChunks = (int) Math.ceil((1.0 * seq.length()) / maxLength);
	             int chunkSize = (int) Math.ceil((1.0 * seq.length()) / numChunks);
	             int pos = 1;
	             while (pos < seq.length()) {
	                 int maxPos = Math.min(seq.length(), pos + chunkSize - 1);
	                 if ((maxPos - pos) >= minLength) {
	                	  seqs.add(new SimpleSequence(
		                             seq.subList(pos, maxPos),
		                             null,
		                             seq.getName() + "_" + pos + "-" + maxPos,
		                             Annotation.EMPTY_ANNOTATION));
	                 }
	                 pos = maxPos + 1;
	             }
	         } else if (seq.length() > minLength) {
	        	 RichSequence.IOTools.writeFasta(System.out, seq, null);
	         }
	     }
	     
	     return seqs;
	}
}
