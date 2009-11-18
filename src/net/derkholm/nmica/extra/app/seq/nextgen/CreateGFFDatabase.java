package net.derkholm.nmica.extra.app.seq.nextgen;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;

import org.biojava.bio.program.gff.GFFRecord;
import org.bjv2.util.cli.Option;

public class CreateGFFDatabase {
	private String dbHost;
	private String dbUser;
	private String dbPassword;
	private String database;
	private boolean createDatabase;
	private boolean dropDatabase;
	private Connection connection;
	private String featureName;
	
	@Option(help="Database host")
	public void setHost(String str) {
		this.dbHost = str;
	}

	@Option(help="Database username")
	public void setUser(String str) {
		this.dbUser = str;
	}

	@Option(help="Database password")
	public void setPassword(String str) {
		this.dbPassword = str;
	}

	@Option(help="Database schema name")
	public void setDatabase(String str) {
		this.database = str;
	}
	
	@Option(help="Feature name")
	public void setFeatureName(String str) {
		this.featureName = str;
	}
	
	

	@Option(help="Create the database, not just the tables", optional=true)
	public void setCreateDatabase(boolean b) {
		this.createDatabase = b;
	}

	@Option(help="Drop database", optional=true)
	public void setDropDatabase(boolean b) {
		this.dropDatabase = b;
	}
	
	private Connection connection() throws Exception {
		if (this.connection == null) {
			if (this.createDatabase) {
				this.connection =
					CountDepths.mysqlConnection(
						dbHost,
						"",
						dbUser,
						dbPassword);

				PreparedStatement statement;
				if (dropDatabase) {
					statement = this.connection.prepareStatement("DROP DATABASE IF EXISTS " + this.database + ";");
					statement.executeUpdate();
					statement.close();
				}
				statement = this.connection.prepareStatement("CREATE DATABASE IF NOT EXISTS " + this.database + ";");
				statement.executeUpdate();
				statement.close();
				statement = this.connection.prepareStatement("USE " + this.database + ";");
				statement.executeUpdate();
				statement.close();
			} else {
				this.connection =
					CountDepths.mysqlConnection(
						dbHost,
						database,
						dbUser,
						dbPassword);

			}
		}
		return this.connection;
	}
	
	public void main(String[] args) throws SQLException, Exception {
		createGFFTable(this.connection(),featureName);
		
	}
	
	private PreparedStatement insertGFFRecordStatement() throws SQLException, Exception {
		/*
		contig_id    	The name of the sequence to which a feature is attached (may actually be a contig, clone, or chromosome name).
		start    		The minimum sequence position covered by the feature
		end    			The maximum position covered by the feature
		strand    		The strand of the feature (should be -1, 0, or 1).
		id    			A unique ID for each feature
		gff_feature    	The "type" of the feature
		gff_source    	The "source" of the feature (e.g. the name of the program which performed the analysis)
		*/
		PreparedStatement insertStatement = 
			connection().prepareStatement(
					String.format("INSERT INTO %s (contig_id,start,end,strand,id,gff_feature,gff_source) " +
					"VALUES (?,?,?,?,?,?)", featureName));
		
		return insertStatement;
	}
	
	public void createGFFTable(Connection conn, String featureName) throws SQLException {
		Statement stat = conn.createStatement();
		stat.executeUpdate(String.format("DROP TABLE if exists %s;",featureName));
		stat.executeUpdate(
				String.format("CREATE TABLE %s (", featureName)
			  +"contig_id    varchar(40) NOT NULL default '',"
			  +"start        int(10) NOT NULL default '0',"
			  +"end          int(10) NOT NULL default '0',"
			  +"strand       int(2) NOT NULL default '0',"
			  +"id           varchar(40) NOT NULL default ''," 
			  +"score        double(16,4) NOT NULL default '0.0000',"
			  +"gff_feature  varchar(40) default NULL,"
			  +"gff_source   varchar(40) default NULL,"
			  +"name         varchar(40) default NULL,"
			  +"hstart       int(11) NOT NULL default '0',"
			  +"hend         int(11) NOT NULL default '0',"
			  +"hid          varchar(40) NOT NULL default'',"
			  +"evalue       varchar(40) default NULL,"
			  +"perc_id      int(10) default NULL,"
			  +"phase        int(11) NOT NULL default '0',"
			  +"end_phase    int(11) NOT NULL default '0',"
			  +"KEY id_contig(contig_id),"
			  +"KEY id_pos(id,start,end)"
			  +") ENGINE=InnoDB;");
		stat.close();
	}
}