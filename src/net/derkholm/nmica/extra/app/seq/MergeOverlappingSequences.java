/*
 * Created on Jul 6, 2005
 */
package net.derkholm.nmica.extra.app.seq;

import gfftools.GFFUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@App(overview = "Write out sequence regions covered by GFF features", 
		generateStub = true)
@NMExtraApp(launchName = "nmmergeoverlapseq", vm = VirtualMachine.SERVER)
public class MergeOverlappingSequences {
	private File featuresFile;
	private Format format = Format.GFF;
	private File seqFile;
	private HashSequenceDB seqDB;
	
	private static enum Format {
		GFF,
		FASTA
	};
	
	@Option(help="Output format")
	public void setFormat(Format format) {
		this.format = format;
	}
	
	@Option(help="Input sequences (only used with fasta output format)", optional=true)
	public void setSeqs(File f) throws FileNotFoundException {
		this.seqFile = f;
		this.seqDB = new HashSequenceDB();
		
		RichSequence.IOTools.readFastaDNA(
			new BufferedReader(new FileReader(f)), null);
	}
	
	@Option(help="The features that are filtered and output " +
			"(those that overlap the masking features are output)")
	public void setFeatures(File f) {
		this.featuresFile = f;
	}
		
    public void main(String[] args) throws Exception {
        final Map<String,Location> locsByChr = GFFUtils.gffToLocationMap(featuresFile);
        
		GFFWriter gffw = null;
		
		if (format == Format.GFF) {
			gffw = new GFFWriter(new PrintWriter(new OutputStreamWriter(System.out)));
		}

        for (String chromo : locsByChr.keySet()) {
        	Location loc = locsByChr.get(chromo);
        	Iterator it = loc.blockIterator();
        	
        	while (it.hasNext()) {
        		Location l = (Location) it.next();
        		
        		if (gffw != null) {
					org.biojava.bio.seq.StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
	        		GFFRecord rec = new SimpleGFFRecord(
							chromo,
							"nmmergeoverlapseq",
							"feature",
							l.getMin(),
							l.getMax(),
							Double.NaN,
							strand,
							0,
							null,
							new HashMap<Object, Object>());
	        		gffw.recordLine(rec);
	        		gffw.endDocument();
        		} else {
        			SymbolList sl = seqDB.getSequence(chromo).subList(l.getMin(), l.getMax());
        			
        		}
        	}
        }
        if (gffw != null) {
            gffw.endDocument();        	
        }
    }
}
