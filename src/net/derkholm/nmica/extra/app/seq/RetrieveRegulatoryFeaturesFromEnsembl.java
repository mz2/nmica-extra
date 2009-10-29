package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.sql.DataSource;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.bio.Annotation;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.JDBCPooledDataSource;
import org.biojavax.bio.seq.RichSequence;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Get regulatory features from Ensembl", generateStub = true)
@NMExtraApp(launchName = "nmensemblregfeat", vm = VirtualMachine.SERVER)
public class RetrieveRegulatoryFeaturesFromEnsembl extends RetrieveEnsemblSequences {

	protected enum Format {
		GFF,
		FASTA
	}
	
	protected Format outputFormat = Format.GFF;
	private String outFileName;
	
	@Option(help="Output format (either gff or fasta, default=gff)",optional=true)
	public void setFormat(Format format) {
		this.outputFormat = format;
	}
	
	@Option(help="Output file (outputs to stdout if not specified)",optional=true)
	public void setOut(String str) {
		this.outFileName = str;
	}
		
	public String schemaBuild() {
		Matcher matcher = Pattern.compile("(\\d+\\_[\\w,\\d]+)$").matcher(this.database);
		
		matcher.find();
		String str = matcher.group(1);
		if (str == null) {
			throw new IllegalArgumentException("The database name doesn't contain a schema build");
		}
		
		System.err.printf("schema build: %s%n",str);
		return str;
	}
	
	private java.sql.PreparedStatement regulatoryFeaturesStatement;
	private java.sql.Connection funcGenConnection;

	protected PreparedStatement regulatoryFeaturesStatement() throws SQLException, Exception {
		if (this.regulatoryFeaturesStatement == null) {
			String str = "SELECT seq_region.name,feature_set.name,seq_region_start,seq_region_end,seq_region_strand FROM regulatory_feature " +
			"LEFT JOIN feature_set ON feature_set.feature_set_id=regulatory_feature.feature_set_id " +
			"LEFT JOIN seq_region ON seq_region.seq_region_id=regulatory_feature.seq_region_id " +
			"WHERE schema_build=?";
			System.err.println(str);
			System.err.println(this.schemaBuild());
			this.regulatoryFeaturesStatement = this.funcGenConnection.prepareStatement(str);
		} 
		
		return this.regulatoryFeaturesStatement;
	}
	
	public void main(String[] args) throws SQLException, Exception {
		String dbURL = String.format("jdbc:mysql://%s:%d/%s", 
											this.host,
											this.port, 
											this.database.replace("core", "funcgen"));

		DataSource db = JDBCPooledDataSource.getDataSource(
				"org.gjt.mm.mysql.Driver", dbURL, username, password);
		this.funcGenConnection = db.getConnection();
		
		regulatoryFeaturesStatement().setString(1, this.schemaBuild());
		regulatoryFeaturesStatement().execute();
		ResultSet resultSet = regulatoryFeaturesStatement().getResultSet();
		
		GFFWriter writer;
		OutputStream os;
		if (this.outputFormat.equals(Format.GFF)) {
			os = null;
			if (this.outFileName == null) {
				writer = new GFFWriter(new PrintWriter(System.out));
			} else {
				writer = new GFFWriter(new PrintWriter(new FileWriter(new File(this.outFileName))));
			}			
		} else {
			writer = null;
			
			if (this.outFileName == null) {
				os = System.out;
			} else {
				os = new BufferedOutputStream(new FileOutputStream(new File(this.outFileName)));
			}
			
			initializeEnsemblConnection();
		}
		
		List<GFFRecord> records = new ArrayList<GFFRecord>();
		while (resultSet.next()) {
			if (this.outFileName != null) System.err.printf("."); //progress indication shown if not outputting to stdout
			String seqRegName = resultSet.getString(1);
			String featSetName = resultSet.getString(2);
			int startCoord = resultSet.getInt(3);
			int endCoord = resultSet.getInt(4);
			//StrandedFeature.Strand strand = (resultSet.getInt(5) == 1) ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE;
			
			if (this.outputFormat.equals(Format.GFF)) {
				SimpleGFFRecord rec = new SimpleGFFRecord();
				rec.setSeqName(seqRegName);
				rec.setFeature(featSetName);
				rec.setStart(startCoord);
				rec.setEnd(endCoord);
				rec.setSource("nmensemblregfeat");
				rec.setStrand(StrandedFeature.UNKNOWN);
				
				writer.recordLine(rec);
				writer.endDocument();
			} else {
				Sequence seq = this.seqDB.getSequence(seqRegName);
				SymbolList sublist = seq.subList(startCoord-1, endCoord-1);
				
				Sequence s = new SimpleSequence(sublist,
													null,
													String.format("%s;%d-%d",
																	seqRegName,
																	startCoord,
																	endCoord), 
																	Annotation.EMPTY_ANNOTATION);
				
				RichSequence.IOTools.writeFasta(os, s, null);
				os.flush();
			}
		}
		writer.endDocument();
	}
}
