package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.extra.app.seq.nextgen.CountDepths.Format;

import org.biojava.bio.BioError;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@NMExtraApp(launchName = "ngmkdepthdb", vm = VirtualMachine.SERVER)
@App(overview = "Create a database to hold sequencing depth pileup data", generateStub = true)
public class CreateDepthDatabase {
	
	private String dbHost;
	private String dbUser;
	private String dbPassword;
	private String database;
	private Connection connection;
	private Format format;
	private File outputFile;
	private boolean createDatabase;

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
	
	@Option(help="Create the database, not just the tables", optional=true)
	public void setCreateDatabase(boolean b) {
		this.createDatabase = b;
	}
	
	private Connection connection() throws Exception {
		if (this.connection == null) {
			if (this.format == Format.MYSQL) {
				this.connection = CountDepths.mysqlConnection(dbHost,database,dbUser,dbPassword);
			} else {
				this.connection = CountDepths.connection(this.format, this.outputFile);
			}
		}
		return this.connection;
	}
	
	public void main(String[] args) throws SQLException, Exception {
		if (this.createDatabase) {
			PreparedStatement statement = connection().prepareStatement("CREATE DATABASE " + this.database);
			statement.executeUpdate();
		}
		
		System.err.println("Creating output tables...");
		if ((format == Format.HSQLDB) || (format == Format.SQLITE) || format == Format.MYSQL) {
			CountDepths.createDepthTable(this.connection());
			CountDepths.createRefSeqTable(this.connection());
		} else {
			throw new BioError("Unsupported format for database creation");
		}
	}

}
