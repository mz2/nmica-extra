package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.utils.ParserException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Remove duplicate entries. As a side effect will also sort the input according to starting position", generateStub=true)
@NMExtraApp(launchName="nmuniqregions", vm=VirtualMachine.SERVER)
public class UniqueRegions {

	private boolean sort;
	private File featuresFile;
	private boolean ignoreStrand;

	@Option(help="Input features (read from stdin if no argument given)",optional=true)
	public void setFeatures(File f) {
		this.featuresFile = f;
	}
	
	@Option(help="Sort the input",optional=true)
	public void setSort(boolean b) {
		this.sort = b;
	}
	
	@Option(help="Ignore strand (default=true)", optional=true)
	public void setIgnoreStrand(boolean b) {
		this.ignoreStrand = b;
	}
	
	public void main(String[] args) throws IOException, BioException, ParserException {
		GFFParser parser = new GFFParser();
		BufferedReader reader;
		
		if (this.featuresFile == null) {
			reader = new BufferedReader(new InputStreamReader(System.in));
		} else {
			reader = new BufferedReader(new FileReader(this.featuresFile));
		}
		
		final boolean ignoreStrand = this.ignoreStrand;
		Comparator<GFFRecord> startPosComparator = new Comparator<GFFRecord>() {

			public int compare(GFFRecord o1, GFFRecord o2) {
				int nameComp = o1.getSeqName().compareTo(o2.getSeqName());
				
				if (nameComp != 0) return nameComp;
				
				int startComp = new Integer(o1.getStart()).compareTo(o2.getStart());
				if (startComp != 0) return startComp;
				
				int endComp = new Integer(o1.getEnd()).compareTo(o2.getEnd());
				if (endComp != 0) return endComp;
				
				if (ignoreStrand) {
					return 0;						
				} else {
					return new Integer(o1.getStrand().getValue()).compareTo(o2.getStrand().getValue());
				}
			}
			
		};
		
		final SortedSet<GFFRecord> recs = new TreeSet<GFFRecord>(startPosComparator);
	
		parser.parse(reader, new GFFDocumentHandler(){
			public void commentLine(String arg0) {
				
			}

			public void endDocument() {
				
			}

			public void recordLine(GFFRecord rec) {
				recs.add(rec);
			}

			public void startDocument(String arg0) {
				
			}
		});
		
		GFFWriter writer = new GFFWriter(new PrintWriter(System.out));
		for (GFFRecord rec : recs) {
			writer.recordLine(rec);
			writer.endDocument();
		}
		writer.endDocument();
	}
}
