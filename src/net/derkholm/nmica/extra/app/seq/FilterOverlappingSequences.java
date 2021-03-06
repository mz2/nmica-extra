/*
 * Created on Jul 6, 2005
 */
package net.derkholm.nmica.extra.app.seq;

import gfftools.GFFUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Map;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

import biobits.utils.IOTools;

@App(overview = "Write out sequence regions covered by GFF features", generateStub = true)
@NMExtraApp(launchName = "nmoverlapseq", vm = VirtualMachine.SERVER)
public class FilterOverlappingSequences {
	private boolean ignoreNames = false;
	//private boolean fraction = false;
	private File featuresFile;
	private File maskFile;
	private Format format = Format.GFF;
	//private boolean mergeOverlapping;
	private boolean exclude;
	private boolean validate;
	
	private static enum Format {
		GFF,
		FRACTION
	};
	
	@Option(help="The features that are filtered and output " +
			"(those that overlap the masking features are output)")
	public void setFeatures(File f) {
		this.featuresFile = f;
	}
	
	@Option(help="The features that are used as a mask")
	public void setMask(File f) {
		this.maskFile = f;
	}
	
	@Option(help="Exclude those sequence regions that overlap with the mask", optional=true)
	public void setExclude(boolean b) {
		this.exclude = b;
	}
	/*
	@Option(help="Merge overlapping features")
	public void setMergeOverlapping(boolean b) {
		this.mergeOverlapping = b;
	}*/
	
	@Option(help="Output the fraction of features covered by the mask",optional=true)
	public void setFraction(boolean b) {
		if (b) {
			this.format = Format.FRACTION;			
		} else {
			this.format = Format.GFF;
		}
	}
	
	@Option(help="Ignore feature names",userLevel=UserLevel.EXPERT,optional=true)
	public void setIgnoreNames(boolean b) {
		this.ignoreNames = b;
	}
	
	@Option(help="Validate input (check that sequence identifiers match)", optional=true)
	public void setValidate(boolean b) {
		this.validate = b;
	}
	
    public void main(String[] args) throws Exception {
        final Map<String,Location> locsByChr = GFFUtils.gffToLocationMap(maskFile);
        final Map<String,Location> featureLocsByChr = GFFUtils.gffToLocationMap(featuresFile);
        
        if (!validate) {
        	System.err.println(
        			"Will not validate the correspondence of sequence identifiers " +
        			"between the features file and the mask file.\n" +
        			"You can turn validation on by adding a flag -validate");
        }
        if (validate &! featureLocsByChr.keySet().equals(locsByChr.keySet())) {
        	System.err.println("The sequence region names are not equal. Will exit.");
        	System.err.printf("Features:%s%n",featureLocsByChr.keySet());
        	System.err.printf("Mask:%s%n",locsByChr.keySet());
        	System.exit(1);
        }
        
        final Location allloc;
        if (ignoreNames) {
        	allloc = LocationTools.union(locsByChr.values());
        } else {
        	allloc = null;
        }
        
        GFFParser gffp = new GFFParser();
        BufferedReader gff = IOTools.fileBufferedReader(featuresFile);
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(System.out));
        final GFFWriter gffw = new GFFWriter(pw);
        gffp.parse(gff, new GFFDocumentHandler() {
        	int all = 0, accepted= 0;
        	
            public void recordLine(GFFRecord record) {
                Location refLoc = allloc == null ? locsByChr.get(record.getSeqName()) : allloc;
                if (refLoc == null) {
                	String n = record.getSeqName();
                	if (n.startsWith("chr")) {
                		n = n.substring(3);
                	}
                	refLoc = locsByChr.get(n);
                	if (refLoc == null) {
                		return;
                	}
                }
                int min = record.getStart();
                int max = record.getEnd();
                if (max < min) {
                	min = record.getEnd();
                	max = record.getStart();
                }
                Location loc = new RangeLocation(min, max);
                ++all;
                
                boolean overlaps = LocationTools.overlaps(loc, refLoc);
                if ((!exclude && overlaps) || (exclude &! overlaps)) {
                	if (format == Format.GFF)
                		gffw.recordLine(record);
                }
            }
            
            public void startDocument(String loc) {
                gffw.startDocument(loc);
            }
            
            public void endDocument() {
                gffw.endDocument();
                
                if (format == Format.FRACTION) {
                	System.out.println((1.0 * accepted) / all);
                }
            }
            
            public void commentLine(String comment) {
                gffw.commentLine(comment);
            }
        }, "");
        pw.flush();
    }
}
