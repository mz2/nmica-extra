package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

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
	
	private Connection connection() throws SQLException {
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
	
	public void main(String[] args) throws SQLException {
		CountDepths.createDepthDatabase(connection());
		
		for (String inFileName : args) {
			
		}
	}
}
