package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.StringTokenizer;

import javax.sql.DataSource;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.biojava.utils.JDBCPooledDataSource;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@App(overview="Follow XREF links to try and find a canonical gene ID", generateStub=true)
@NMExtraApp(launchName = "nmcanonicalizegenes", vm = VirtualMachine.SERVER)
public class CanonicalizeGeneList {

	/*
	private String dbURL;
	private String dbUser;
	private String dbPass = "";
	*/
	private boolean useDisplayLabels = false;
	private String username = "anonymous";
	private String password = "";
	private String host = "ensembldb.ensembl.org";
	private String database;
	private int port = 5306;
	private int schemaVersion = 54;
	
	private PreparedStatement get_gsiForTranslation;
	private Connection connection;
	private PreparedStatement get_geneStableId;
	private PreparedStatement get_xref;
	private PreparedStatement get_gsiForTranscript;
	
	@Option(help="...", optional=true)
	public void setUseDisplayLabels(boolean  b) {
		this.useDisplayLabels = b;
	}
	
	@Option(help = "Ensembl password (default='')", optional = true)
	public void setPassword(String password) {
		this.password = password;
	}

	@Option(help = "Ensembl hostname (default=ensembldb.ensembl.org", optional = true)
	public void setHost(String host) {
		this.host = host;
	}
	
	@Option(help = "Ensembl database to retrieve sequences from (e.g. 'danio_rerio_core_54_8')")
	public void setDatabase(String db) {
		this.database = db;
	}

	@Option(help = "Ensembl username (default=anonymous)", optional = true)
	public void setUser(String user) {
		this.username = user;
	}

	public Connection connection() throws Exception {
		if (this.connection == null) {
			DataSource db = JDBCPooledDataSource.getDataSource(
	                "org.gjt.mm.mysql.Driver",
	                String.format("jdbc:mysql://%s:%d/%s", this.host,
	        				this.port, this.database),
	                username,
	                password);
			this.connection = db.getConnection();
		}
		return connection;
	}

	private void initPreparedStatements() throws SQLException, Exception {
		this.get_geneStableId = connection().prepareStatement(
				"select gene_id from gene_stable_id where stable_id = ?"
	    );
		
		if (useDisplayLabels) {
			this.get_xref = connection().prepareStatement(
				"select ensembl_id, ensembl_object_type " +
				"  from xref, object_xref " +
				" where object_xref.xref_id = xref.xref_id and " +
				"       xref.display_label = ?"
			);
		} else {
			this.get_xref = connection().prepareStatement(
					"select ensembl_id, ensembl_object_type " +
					"  from xref, object_xref " +
					" where object_xref.xref_id = xref.xref_id and " +
					"       xref.dbprimary_acc = ?"
		    );
		}
		this.get_gsiForTranscript = connection().prepareStatement(
				"select stable_id " +
				"  from gene_stable_id, transcript " +
				" where gene_stable_id.gene_id = transcript.gene_id and " +
				"       transcript.transcript_id = ?"
		);
		this.get_gsiForTranslation = connection().prepareStatement(
				"select stable_id " +
				"  from gene_stable_id, transcript, translation " +
				" where gene_stable_id.gene_id = transcript.gene_id and " +
				"       transcript.transcript_id = translation.transcript_id and " +
				"       translation.translation_id = ?"
		);


	}
	/**
	 * @param args
	 */
	public void main(String[] args) 
		throws Exception
	{
		initPreparedStatements();
		
		Reader r;
		if (args.length > 0) {
			r = new FileReader(args[0]);
		} else {
			r = new InputStreamReader(System.in);
		}
		
		BufferedReader br = new BufferedReader(r);
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			line = new StringTokenizer(line.trim(), " ,").nextToken();
			
			{
				get_geneStableId.setString(1, line);
				ResultSet rs = get_geneStableId.executeQuery();
				boolean success = rs.next();
				rs.close();
				if (success) {
					System.out.println(line);
					continue;
				}
			}
			{
				get_xref.setString(1, line);
				ResultSet rs = get_xref.executeQuery();
				int tid = -1;
				int tlid = -1;
				while (rs.next()) {
					int id = rs.getInt(1);
					String type = rs.getString(2);
					if (type.equalsIgnoreCase("transcript")) {
						tid = id;
					} else if (type.equalsIgnoreCase("translation")) {
						tlid = id;
					}
				}
				rs.close();
				
				if (tid > 0) {
					get_gsiForTranscript.setInt(1, tid);
					rs = get_gsiForTranscript.executeQuery();
					rs.next();
					System.out.println(rs.getString(1) + "\t" + line);
					rs.close();
				} else if (tlid > 0) {
					get_gsiForTranslation.setInt(1, tlid);
					rs = get_gsiForTranslation.executeQuery();
					rs.next();
					System.out.println(rs.getString(1) + "\t" + line);
					rs.close();
				} else {
					// System.out.println("Skipping " + line);
				}
			}
		}
	}

}
