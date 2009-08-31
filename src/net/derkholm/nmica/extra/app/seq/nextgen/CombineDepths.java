package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@NMExtraApp(launchName = "ngcombinedepths", vm = VirtualMachine.SERVER)
@App(overview = "Combine depths to a single SQLite database", generateStub = true)
public class CombineDepths {
	private File outputFile;
	private Connection connection;

	@Option(help="Output file")
	public void setOut(File f) {
		this.outputFile = f;
	}
	
	private Connection connection() throws SQLException, ClassNotFoundException {
		Class.forName("org.sqlite.JDBC");
		if (this.connection == null) {
			this.connection = 
				DriverManager.getConnection(
					String.format(
						"jdbc:sqlite:%s",
						this.outputFile.getPath()));

			this.connection.setAutoCommit(false);
		}
		return this.connection;
	}
	
	public void main(String[] args) throws SQLException, ClassNotFoundException {
		CountDepths.createDepthDatabase(connection());
		
		for (String inFileName : args) {
			System.err.printf("Adding %s...%n", inFileName);
			Connection conn = 
				DriverManager.getConnection(String.format(
					"jdbc:sqlite:%s",inFileName));
			
			Statement statement = conn.createStatement();
			ResultSet results = statement.executeQuery(
				"SELECT ref_name,begin_coord,end_coord,depth,pvalue FROM window;");
			
			while (results.next()) {
				String refName = results.getString(1);
				int begin = results.getInt(2);
				int end = results.getInt(3);
				double d = results.getDouble(4);
				double p = results.getDouble(5);
				
				System.err.printf("%s %d %d %.2f %.2f%n",refName,begin,end,d,p);
			}
			results.close();
		}
	}
}
