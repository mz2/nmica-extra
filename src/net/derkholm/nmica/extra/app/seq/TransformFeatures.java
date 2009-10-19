package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.sql.SQLException;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.bjv2.util.cli.App;

@App(overview="Transform features in GFF files (move / scale)", generateStub=true)
@NMExtraApp(launchName="nmtransformfeat", vm=VirtualMachine.SERVER)
public class TransformFeatures {

	private File featuresFile;
	private int scaleTo;
	private int moveBy;
	private String outFile;
	private GFFWriter gffw;

	public void setFeatures(File f) {
		this.featuresFile = f;
	}
	
	public void setScaleTo(int i) {
		this.scaleTo = i;
	}
	
	public void setMoveBy(int i) {
		this.moveBy = i;
	}
	
	public void setOut(String f) {
		this.outFile = f;
	}
	
	public void main(String[] args) throws SQLException, Exception {
		
		final OutputStream os;
		if (this.outFile == null) {
			os = System.out;
		} else {
			os = new FileOutputStream(this.outFile);
		}

		InputStream inputStream;
		if (featuresFile == null) {
			inputStream = System.in;
		} else {
			inputStream = new FileInputStream(this.featuresFile);
		}

		
		gffw = new GFFWriter(new PrintWriter(new OutputStreamWriter(os)));
		
		GFFParser parser = new GFFParser();
		parser.parse(
				new BufferedReader(new InputStreamReader(inputStream)),
				new GFFDocumentHandler() {

					@Override
					public void commentLine(String arg0) {
						
					}

					@Override
					public void endDocument() {
						
					}

					@Override
					public void recordLine(GFFRecord r) {
						System.err.printf(".");

						
						GFFRecord rec = new SimpleGFFRecord(
								r.getSeqName(),
								r.getSource(),
								r.getFeature(),
								r.getStart() + moveBy,
								r.getEnd() + moveBy,
								r.getScore(),
								r.getStrand(),
								r.getFrame(),
								r.getComment(),
								r.getGroupAttributes());

						gffw.recordLine(rec);
						gffw.endDocument();
						
					}

					@Override
					public void startDocument(String arg0) {
						// TODO Auto-generated method stub
						
					}
					
					
				});
	}	
}