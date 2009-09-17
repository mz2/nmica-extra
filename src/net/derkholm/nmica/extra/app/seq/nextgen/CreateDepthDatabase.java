package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.extra.app.seq.nextgen.CountDepths.Format;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
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
	private Format format = Format.MYSQL;
	private File outputFile;
	private boolean createDatabase = false;
	private boolean dropDatabase = false;
	private boolean dropTable = false;
	private List<String> refSeqNames;
	private Map<String, Integer> refLengths;
	private HashMap<String, Integer> refIds;
	private Map<String, Integer> readCounts;

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

	@Option(help="Output format (default=mysql)", optional=true)
	public void setFormat(Format format) {
		this.format = format;
	}

	@Option(help="Create the database, not just the tables", optional=true)
	public void setCreateDatabase(boolean b) {
		this.createDatabase = b;
	}

	@Option(help="Drop and recreate database", optional=true)
	public void setDropDatabase(boolean b) {
		this.dropDatabase = b;
	}
	
	@Option(help="Drop tables before inserting", optional=true)
	public void setDropTable(boolean b) {
		this.dropTable = b;
	}


	@Option(help = "Reference sequence lengths")
	public void setRefLengths(File f) throws BioException, IOException {
		try {
			this.refLengths = SAMProcessor.parseRefLengths(f);
		} catch (Exception e) {
			e.printStackTrace();
		}
		this.refSeqNames = SAMProcessor.parseRefNamesFromRefLengthFile(f);
	}

	@Option(help = "Read counts")
	public void setReadCounts(File f) {
		try {
			this.readCounts = SAMProcessor.parseReadCounts(f);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private Connection connection() throws Exception {
		if (this.connection == null) {
			if (this.format == Format.MYSQL) {
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
			} else {
				this.connection =
					CountDepths.connection(
						this.format,
						this.outputFile);
			}
		}
		return this.connection;
	}

	public void main(String[] args) throws SQLException, Exception {
		System.err.println("Creating output tables...");
		if ((format == Format.HSQLDB) || (format == Format.SQLITE) || format == Format.MYSQL) {
			CountDepths.createDepthTable(this.connection(),this.dropTable);
			CountDepths.createRefSeqTable(this.connection(),this.dropTable);
			CreateDepthDatabase.createConservationTable(this.connection(),this.dropTable);
		} else {
			throw new BioError("Unsupported format for database creation");
		}

		this.refIds = new HashMap<String, Integer>();

		int i = 0;

		PreparedStatement stat;
		try {
			stat = CountDepths.insertRefSeqNameStatement(this.connection());

			for (String name : this.refSeqNames) {
				this.refIds.put(name, i);
				stat.setInt(1, i);
				stat.setString(2, name);
				stat.setDouble(3, (double)this.readCounts.get(name));
				stat.addBatch();
				stat.executeBatch();
				stat.clearBatch();
				connection().commit();
				i++;
			}
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public static void createConservationTable(Connection conn, boolean dropTable) throws SQLException {
		Statement stat = conn.createStatement();

		if (dropTable) {
			stat.executeUpdate("DROP TABLE IF EXISTS conservation;");
		}
		stat.executeUpdate("CREATE TABLE IF NOT EXISTS conservation ("
							+ "id integer,"
							+ "ref_id INTEGER,"
							+ "coord INTEGER,"
							+ "conservation DOUBLE,"
							+ " PRIMARY KEY (id));");
		stat.executeUpdate("CREATE INDEX ref_name_begin_idx ON conservation(ref_id,coord);");
		stat.close();
	}
}
