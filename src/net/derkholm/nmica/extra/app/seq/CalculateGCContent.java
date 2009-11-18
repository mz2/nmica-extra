package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.extra.app.seq.nextgen.SAMProcessor;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@NMExtraApp(launchName = "nmgccontent", vm = VirtualMachine.SERVER)
@App(overview = "Calculate the GC content of sequence windows", generateStub = true)
public class CalculateGCContent {

	private File seqsFilename;
	private int winSize;
	private File out;

	@Option(help="Sequence filename")
	public void setSeqs(File f) {
		this.seqsFilename = f;

	}
	
	@Option(help="Window size")
	public void setWindowSize(int winSize) {
		this.winSize = winSize;
	}
	
	@Option(help="Output file")
	public void setOut(File f) {
		this.out = f;
	}
	
	private enum Format {
		TSV,
		SQLITE
	};
	
	private Format format = Format.SQLITE;
	private Connection connection;
	private HashMap<String, Integer> refIds;
	private Map<String, Integer> refSeqLengths;
	private List<String> refSeqNames;
	
	public void setFormat(Format format) {
		this.format = format;
	}
	
	@Option(help="Reference sequence lengths")
	public void setRefSeqLengths(File f) {
		this.refSeqLengths = SAMProcessor.parseRefLengths(f);
		this.refSeqNames = SAMProcessor.parseRefNamesFromRefLengthFile(f);
	}
	
	public static void createGCContentTable(Connection conn) throws SQLException {
		Statement stat = conn.createStatement();
		stat.executeUpdate("DROP TABLE if exists gccontent;");
		stat.executeUpdate(
			"CREATE TABLE gccontent (" +
				"id integer primary key autoincrement," +
				"ref_id integer," +
				"begin_coord integer," +
				"end_coord integer," +
				"gccontent float) ENGINE=InnoDB;");
		stat.executeUpdate("CREATE INDEX gc_begin_end ON gccontent(ref_id,begin_coord,end_coord);");
		stat.close();
	}
	
	private Connection connection() throws SQLException {
		if (this.connection == null) {
			this.connection = 
				DriverManager.getConnection(
					String.format(
						"jdbc:sqlite:%s",
						this.out.getPath()));

			this.connection.setAutoCommit(true);
		}
		return this.connection;
	}
	
	public void main(String[] args) throws FileNotFoundException, BioException, SQLException, ClassNotFoundException {
		Class.forName("org.sqlite.JDBC");
		createGCContentTable(connection());
		BufferedReader br = new BufferedReader(new FileReader(seqsFilename));
		SequenceIterator stream = SeqIOTools.readFastaDNA(br);
		
		PreparedStatement insertStatement = connection().prepareStatement("INSERT INTO gccontent VALUES (?, ?, ?, ?, ?)");
		int id = 0;
		
		while (stream.hasNext()) {
			Sequence seq = stream.nextSequence();
			System.err.printf("Calculating GC content for %s%n");

			int refId = getRefId(seq.getName());
			
			int len = seq.length();
			
			for (int i = 1; i <= (len-winSize); i++) {
				SymbolList symList = seq.subList(i, i+winSize);
				
				/* ambiguous positions are ignored */
				int gc = 0;
				int at = 0;
				
			    for (int pos = 1; pos <= winSize; ++pos) {
					Symbol sym = symList.symbolAt(pos);
					if (sym.equals(DNATools.g()) || sym.equals(DNATools.c())) {++gc;}
					else if (sym.equals(DNATools.a()) || sym.equals(DNATools.t())) {++at;}
			    }
			    double gcRatio = (double)gc/(double)(gc+at);
			   
			    /*
			    if (!Double.isNaN(gcRatio)) {
				    insertStatement.setInt(1, id++);
				    insertStatement.setInt(2, refId);
				    insertStatement.setInt(3, i);
				    insertStatement.setInt(4, i+winSize);
				    insertStatement.setFloat(5, (float) gcRatio);
				    insertStatement.addBatch();
				    insertStatement.executeBatch();
			    }*/
			}
			
			break;
		}
	}

	public int getRefId(String seqName) {
		if (this.refIds == null) {
			this.refIds = new HashMap<String,Integer>();
			
			int i = 0;
			for (String name : this.refSeqNames) {
				this.refIds.put(name, i++);
			}
		}
		return refIds.get(seqName);
	}
}