package net.derkholm.nmica.extra.app.seq;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.sql.SQLException;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Get noncoding sequences from Ensembl for motif discovery", generateStub = true)
@NMExtraApp(launchName = "nmensemblchromosomes", vm = VirtualMachine.SERVER)
public class RetrieveChromosomesFromEnsembl extends RetrieveEnsemblSequences {
	
	private File outFile;

	@Option(help="Output file",optional=true)
	public void setOut(File f) {
		this.outFile = f;
	}

	public void main(String[] args) throws SQLException, Exception {
		initializeEnsemblConnection();
		
		PrintStream os = null;
		if (this.outFile == null) {
			os = System.out;
		} else {
			os = new PrintStream(new FileOutputStream(this.outFile, true));
		}
		
		SequenceDB chromos = ensemblConnection.getSequenceDB("chromosome");
		for (Object id : chromos.ids()) {
			String idStr = (String)id;
			
			Sequence seq = chromos.getSequence(idStr);
			System.err.printf("Retrieved chromosome %s (length : %d)%n", idStr, seq.length());
			RichSequence.IOTools.writeFasta(os, seq, null);
		}
		
		os.close();
	}
}
