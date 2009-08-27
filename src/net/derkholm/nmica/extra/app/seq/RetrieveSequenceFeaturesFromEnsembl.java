package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.sql.SQLException;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Get noncoding sequences from Ensembl for motif discovery", generateStub = true)
@NMExtraApp(launchName = "nmensemblfeat", vm = VirtualMachine.SERVER)
public class RetrieveSequenceFeaturesFromEnsembl extends RetrieveEnsemblSequences {
	
	private File outFile;
	private File featuresFile;

	@Option(help="Output file",optional=true)
	public void setOut(File f) {
		this.outFile = f;
	}
	
	@Option(help="Features file (read from stdin if not included)", optional=true) 
	public void setFeatures(File f) {
		this.featuresFile = f;
	}

	public void main(String[] args) throws SQLException, Exception {
		initializeEnsemblConnection();
		
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
		
		GFFParser parser = new GFFParser();
		parser.parse(
				new BufferedReader(new InputStreamReader(inputStream)), 
				new GFFDocumentHandler() {
			public void commentLine(String str) {}
			public void endDocument() {}
			public void startDocument(String str) {}

			public void recordLine(GFFRecord recLine) {
				try {
					SymbolList symList = seqDB.getSequence(recLine.getSeqName()).subList(recLine.getStart(), recLine.getEnd());
					if (recLine.getStrand().equals(StrandedFeature.NEGATIVE)) {
						symList = DNATools.reverseComplement(symList);
					}
					
					Sequence s = new SimpleSequence(symList, null, 
							String.format("%s;%d-%d(%s)", 
									recLine.getSeqName(), 
									recLine.getStart(), 
									recLine.getEnd(), 
									recLine.getStrand()
										.equals(StrandedFeature.POSITIVE)? "+" : "-"),
							Annotation.EMPTY_ANNOTATION);
					
					RichSequence.IOTools.writeFasta(os, s, null);

				} catch (IllegalIDException e) {
					e.printStackTrace();
				} catch (BioException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		});
	}
}