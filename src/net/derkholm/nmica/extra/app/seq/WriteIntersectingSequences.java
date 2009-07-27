package net.derkholm.nmica.extra.app.seq;

import gfftools.GFFUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFTools;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

@App(overview = "Write out sequence regions covered by GFF features", generateStub = true)
@NMExtraApp(launchName = "nmintersectseq", vm = VirtualMachine.SERVER)
public class WriteIntersectingSequences {
    public static enum Format {
    	GFF,
    	FASTA
    }
    
	private File seqsFile;
	private boolean negate;
	private SequenceDB seqDB;
	private net.derkholm.nmica.extra.app.seq.WriteIntersectingSequences.Format outputFormat;
	private boolean validate;

	@Option(help = "Input sequences (output the intersecting sequences rather than a new GFF file)", optional=true)
	public void setSeqs(File f) throws FileNotFoundException, ChangeVetoException, NoSuchElementException, BioException {
		this.seqsFile = f;
		
		this.seqDB = new HashSequenceDB();
		for (SequenceIterator si = RichSequence
				.IOTools
					.readFastaDNA(
						new BufferedReader(
							new FileReader(f)), null); si.hasNext();) {
			
			seqDB.addSequence(si.nextSequence());
		}
	}
	
	@Option(help = "Negate output after intersecting (i.e. substract features rather than intersecting)", optional=true)
	public void setNegate(boolean b) {
		this.negate = b;
	}
	
	@Option(help="Output format: fasta|gff", optional=true)
	public void setFormat(Format format) {
		this.outputFormat = format;
	}
	
	@Option(help="Validate input (check that sequence identifiers match)", optional=true, userLevel = UserLevel.DEBUG)
	public void setValidate(boolean b) {
		this.validate = b;
	}
	
    /**
     * @param args
     * @throws Exception 
     */
    public void main(String[] args) throws Exception {
        List<Map<String,Location>> locs = new ArrayList<Map<String,Location>>();
        for (String fileName : args) {
            locs.add(GFFUtils.gffToLocationMap(new File(fileName)));
        }
        
        Set<String> seqIds;
        {
            Iterator<Map<String,Location>> i = locs.iterator();
            seqIds = new HashSet<String>(i.next().keySet());
            while (i.hasNext()) {
                seqIds.retainAll(i.next().keySet());
            }
        }
        
        if (validate && (seqDB != null)) {
        	for (Map<String, Location> ls : locs) {
            	WriteCoveredSequences.validateGFFSequenceIdentifiersAgainstSequences(ls, seqDB);        		
        	}
        }
        
        PrintWriter pw = null;
        GFFWriter gffw = null;
        GFFEntrySet gffEntries = null;
        
        if (outputFormat == Format.GFF) {
        	pw = new PrintWriter(new OutputStreamWriter(System.out));
        	gffw = new GFFWriter(pw);
        } else {
            gffEntries = new GFFEntrySet();
        }
        
       
        for (String id : seqIds) {
            Iterator<Map<String,Location>> i = locs.iterator();
            Location l = i.next().get(id);
            while (i.hasNext()) {
                l = LocationTools.intersection(l, i.next().get(id));
            }
            
            if (negate) {
				l = LocationTools.subtract(new RangeLocation(1, seqDB.getSequence(id).length()), l);
			}
            
            SimpleGFFRecord r = new SimpleGFFRecord();
            r.setSeqName(id);
            r.setFeature("block");
            r.setSource("nmintersectseq");
            r.setStrand(StrandedFeature.POSITIVE);
            for (Iterator<?> bi = l.blockIterator(); bi.hasNext(); ) {
                Location bloc = (Location) bi.next();
                r.setStart(bloc.getMin());
                r.setEnd(bloc.getMax());
                
                if (gffw != null) {gffw.recordLine(r);}
                else {
                	Sequence seq = new SimpleSequence(seqDB.getSequence(id).subList(bloc.getMin(),bloc.getMax()),null,
                			String.format("%s_|%d-%d|",id,bloc.getMin(),bloc.getMax()),Annotation.EMPTY_ANNOTATION);
                	RichSequence.IOTools.writeFasta(System.out, seq, null);
                }
            }
        }
        if (pw != null) pw.flush();
        
        if (seqDB != null) {
        	System.err.println("Writing output sequences...");
            GFFTools.annotateSequences(seqDB, gffEntries);
        }
    }
}
