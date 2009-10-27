/*
 * Created on Jul 20, 2005
 */
package gfftools;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;

import biobits.utils.Collects;
import biobits.utils.IOTools;

public class GFFUtils {

	public static GFFParser parser() {
		if (Boolean.getBoolean("gff3")) {
			// System.err.println("Making a GFF3 parser");
			return new GFFParser3();
		} else {
			// System.err.println("Making a GFF2 parser");
			return new GFFParser();
		}
	}
	
    private GFFUtils() {
    }
    
    public static Map<String,List<GFFRecord>> gffToRecordMap(File file)
    	throws Exception
    {
    	final Map<String,List<GFFRecord>> blocs = new HashMap<String,List<GFFRecord>>();
        GFFParser gffp = new GFFParser();
        gffp.parse(IOTools.fileBufferedReader(file), new GFFDocumentHandler() {
            public void startDocument(String locator) {
            }
    
            public void endDocument() {
            }
    
            public void commentLine(String comment) {
            }
    
            public void recordLine(GFFRecord record) {
                Collects.pushOntoMap(blocs, record.getSeqName(), record);
            }
        });
        return blocs;
    }
    
    public static Map<String,Location> gffToLocationMap(File file) 
        throws Exception
    {
        final Map<String,List<Location>> blocs = new HashMap<String,List<Location>>();
        GFFParser gffp = new GFFParser();
        gffp.parse(IOTools.fileBufferedReader(file), new GFFDocumentHandler() {
    
            public void startDocument(String locator) {
            }
    
            public void endDocument() {
            }
    
            public void commentLine(String comment) {
            }
    
            public void recordLine(GFFRecord record) {
            	int min = record.getStart();
                int max = record.getEnd();
                if (max < min) {
                	min = record.getEnd();
                	max = record.getStart();
                }
                Location bloc = new RangeLocation(min, max);
                List<Location> locList = blocs.get(record.getSeqName());
                if (locList == null) {
                    locList = new ArrayList<Location>();
                    blocs.put(record.getSeqName(), locList);
                }
                locList.add(bloc);
            }
        });
        
        Map<String,Location> unions = new HashMap<String,Location>();
        for (Map.Entry<String,List<Location>> me : blocs.entrySet()) {
            unions.put(me.getKey(), LocationTools.union(me.getValue()));
        }
        return unions;
    }

	public static String gaAsString(GFFRecord r, String key) {
		return ((List<?>) r.getGroupAttributes().get(key)).get(0).toString();
	}
}
