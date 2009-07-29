/*
 * Created on Jul 6, 2005
 */
package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = 
	"Write out sequence regions covered by GFF features " +
	"in all input files with the specified count of overlapping nucleotides", 
	generateStub = true)
@NMExtraApp(launchName = "nmoverlapseq", vm = VirtualMachine.SERVER)
public class FilterOverlappingSequences {
    private String masterBrandKey = null;
    private String outputBrandKey = null;
    private int minOverlap = 1;
        
    @Option(help="Minimum overlap between features (default=1)")
    public void setMinOverlap(int i) {
        this.minOverlap = i;
    }
    
    public void setMasterBrandKey(String masterBrandKey) {
        this.masterBrandKey = masterBrandKey;
    }

    public void setOutputBrandKey(String outputBrandKey) {
        this.outputBrandKey = outputBrandKey;
    }

    private static final class BrandedSpan implements Comparable<BrandedSpan> {
        public final int min;
        public final int max;
        public final String brand;
        
        private BrandedSpan(int min, int max, String brand) {
            this.min = min;
            this.max = max;
            this.brand = brand;
        }

        public int compareTo(BrandedSpan arg0) {
            int dif = min - arg0.min;
            if (dif != 0) {
                return dif;
            } 
            dif = max - arg0.max;
            if (dif != 0) {
                return dif;
            }
            return brand.compareTo(arg0.brand);
        }        
    }
    
    public void main(String[] args)
        throws Exception
    {
        File ref = new File(args[0]);
        File f = new File(args[1]);
        
        final Map<String,List<BrandedSpan>> locsByChr = new HashMap<String,List<BrandedSpan>>();
        GFFParser gffp = new GFFParser();
        BufferedReader refGff = new BufferedReader(new FileReader(ref));
        gffp.parse(refGff, new GFFDocumentHandler() {
            public void startDocument(String locator) {
            }

            public void endDocument() {
            }

            public void commentLine(String comment) {
            }

            public void recordLine(GFFRecord record) {
                String brand;
                if (masterBrandKey != null) {
                	try {
                		brand = ((List) record.getGroupAttributes().get(masterBrandKey)).get(0).toString();
                	} catch (Exception ex) {return;}
                } else {
                    brand = record.getFeature() + "_" + record.getSeqName() + "_" + record.getStart();
                }
                BrandedSpan b = new BrandedSpan(record.getStart(), record.getEnd(), brand);
                List<BrandedSpan> l = locsByChr.get(record.getSeqName());
                if (l == null) {
                    l = new ArrayList<BrandedSpan>();
                    locsByChr.put(record.getSeqName(), l);
                }
                l.add(b);
            }
        });
        
        for (List<BrandedSpan> l : locsByChr.values()) {
            Collections.sort(l);
        }
        
        
        BufferedReader gff = new BufferedReader(new FileReader(f));
        PrintWriter pw = new PrintWriter(new OutputStreamWriter(System.out));
        final GFFWriter gffw = new GFFWriter(pw);
        gffp.parse(gff, new GFFDocumentHandler() {
            public void recordLine(GFFRecord record) {
                List<BrandedSpan> refLocs = locsByChr.get(record.getSeqName());
                if (refLocs == null) {
                    return;
                }
                
                int lb = 0, ub = refLocs.size() - 1;
                while (ub >= lb) {
                    int mid = (lb + ub) / 2;
                    BrandedSpan midSpan = refLocs.get(mid);
                    if (midSpan.min > record.getEnd()) {
                        ub = mid - 1;
                    } else if (midSpan.max < record.getStart()) {
                        lb = mid + 1;
                    } else {
                        int overlap = Math.min(midSpan.max, record.getEnd()) - Math.max(midSpan.min, record.getStart()) + 1;
                        if (overlap >= minOverlap) {
                            SimpleGFFRecord r = new SimpleGFFRecord(record);
                            r.getGroupAttributes().put(outputBrandKey, Collections.singletonList(midSpan.brand));
                            gffw.recordLine(r);
                        }
                        return;
                    }
                }
            }
            
            public void startDocument(String loc) {
                gffw.startDocument(loc);
            }
            
            public void endDocument() {
                gffw.endDocument();
            }
            
            public void commentLine(String comment) {
                gffw.commentLine(comment);
            }
        }, "");
        pw.flush();
    
    }
}
