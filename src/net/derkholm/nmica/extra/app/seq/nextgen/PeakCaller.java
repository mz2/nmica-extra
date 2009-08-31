package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Map;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.extra.peak.Peak;
import net.derkholm.nmica.extra.peak.Window;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@NMExtraApp(launchName = "ngdepth", vm = VirtualMachine.SERVER)
@App(overview = "Output sequencing depth inside a window.", generateStub = true)
public class PeakCaller {
	
	private Connection depthConnection;
	private Connection controlDepthConnection;
	private File depthFile;
	private File controlDepthFile;
	private Map<String, Integer> readCounts;
	private Map<String, Integer> refSeqLengths;
	private double pvalue = 1e-6;
	private double fdr = 0.1;
	private Connection connection;
	private File outputFile;
	private PreparedStatement depthAtPositionStatement;
	private PreparedStatement controlDepthAtPositionStatement;
	private Peak peak;

	@Option(help="Sequencing depths")
	public void setDepths(File f) throws ClassNotFoundException, SQLException {
		Class.forName("org.sqlite.JDBC");
		this.depthFile = f;
		if (this.depthConnection == null) {
			this.depthConnection = 
				DriverManager.getConnection(
					String.format(
						"jdbc:sqlite:%s",
						f.getPath()));

			this.depthConnection.setAutoCommit(true);
		}
	}
	
	@Option(help="Control sequencing depths")
	public void setControlDepths(File f) throws ClassNotFoundException, SQLException {
		Class.forName("org.sqlite.JDBC");
		this.controlDepthFile = f;
		if (this.controlDepthConnection == null) {
			this.controlDepthConnection = 
				DriverManager.getConnection(
					String.format(
						"jdbc:sqlite:%s",
						f.getPath()));

			this.controlDepthConnection.setAutoCommit(true);
		}
	}
	
	@Option(help="Output file")
	public void setOut(File f) {
		this.outputFile = f;
	}
	/*
	@Option(help="Read counts")
	public void setReadCounts(File f) {
		this.readCounts = SAMProcessor.parseReadCounts(f);
	}*/
	
	@Option(help="Reference sequence lengths")
	public void setRefLengths(File f) {
		this.refSeqLengths = SAMProcessor.parseRefLengths(f);
	}
	
	@Option(help="P-value cutoff (default=1e-6)",optional=true)
	public void setPValue(double d) {
		this.pvalue = d;
	}
	
	@Option(help="FDR cutoff (default=0.10)",optional=true)
	public void setFDR(double d) {
		this.fdr = d;
	}
	
	private void createPeakDatabase() throws SQLException, ClassNotFoundException {
		Statement stat = connection().createStatement();
		stat.executeUpdate("DROP TABLE if exists peak;");
		stat.executeUpdate("CREATE TABLE peak (" +
				"id integer primary key," +
				"ref_name varchar," +
				"begin_coord integer," +
				"end_coord integer," +
				"depth float," +
				"depth_control float," +
				"pvalue float," +
				"fdr float);");
		stat.executeUpdate("CREATE INDEX ref_name_begin_end_idx ON depth(ref_name,begin_coord,end_coord);");
		stat.executeUpdate("CREATE INDEX ref_name_begin_idx ON depth(ref_name,begin_coord);");
	}
	
	private PreparedStatement 
		getDepthAtPositionStatement() throws SQLException {
		if (this.depthAtPositionStatement == null) {
			this.depthAtPositionStatement = 
				depthConnection.prepareStatement(
					"SELECT begin_coord,end_coord,depth,pvalue " +
					"FROM window WHERE ref_name = ? AND begin_coord = ?");
		}
		return this.depthAtPositionStatement;
	}
	
	private PreparedStatement 
		getControlDepthAtPositionStatement() throws SQLException {
		if (this.controlDepthAtPositionStatement == null) {
			this.controlDepthAtPositionStatement = 
				controlDepthConnection.prepareStatement(
					"SELECT begin_coord,end_coord,depth,pvalue " +
					"FROM window WHERE ref_name = ? AND begin_coord = ?");
		}
		return this.controlDepthAtPositionStatement;
	}

	private Window getDepth(String refSeq, int position) throws SQLException {
		this.getDepthAtPositionStatement().setString(1, refSeq);
		this.getDepthAtPositionStatement().setInt(2, position);
		this.getDepthAtPositionStatement().execute();
		ResultSet results = this.getDepthAtPositionStatement().getResultSet();
		int beginCoord = results.getInt(1);
		int endCoord = results.getInt(2);
		double depth = results.getDouble(3);
		double pvalue = results.getDouble(4);
		
		return new Window(refSeq, beginCoord, endCoord, depth, pvalue);
	}
	
	private Window getControlDepth(String refSeq, int position) throws SQLException {
		this.getControlDepthAtPositionStatement().setString(1, refSeq);
		this.getControlDepthAtPositionStatement().setInt(2, position);
		this.getControlDepthAtPositionStatement().execute();
		ResultSet results = this.getControlDepthAtPositionStatement().getResultSet();
		int beginCoord = results.getInt(1);
		int endCoord = results.getInt(2);
		double depth = results.getDouble(3);
		double pvalue = results.getDouble(4);
		
		return new Window(refSeq, beginCoord, endCoord, depth, pvalue);
	}
	
	private Connection connection() throws SQLException, ClassNotFoundException {
		Class.forName("org.sqlite.JDBC");
		if (this.connection == null) {
			this.connection = 
				DriverManager.getConnection(
					String.format(
						"jdbc:sqlite:%s",
						this.outputFile.getPath()));

			this.connection.setAutoCommit(true);
		}
		return this.connection;
	}
	
	public void main(String[] args) throws SQLException, ClassNotFoundException {
		this.createPeakDatabase();
		for (String refSeq : this.refSeqLengths.keySet()) {
			//int readCount = this.readCounts.get(refSeq);
			int refSeqLength = this.refSeqLengths.get(refSeq);
			
			
			for (int i = 0; i < refSeqLength; i++) {
				Window depth = this.getDepth(refSeq, i);
				Window controlDepth = this.getControlDepth(refSeq, i);
				
				if (depth.getPvalue() > this.pvalue) {
					if (this.peakIsOpen()) {
						this.endPeak();
					}
				} else {
					if (!this.peakIsOpen()) {
						this.startPeak(refSeq,i);
					}
				}
			}
		}
	}

	private void startPeak(String refSeq, int i) {
		this.peak = new Peak(refSeq, i);
	}

	private void endPeak() {
		// TODO Auto-generated method stub
		
	}

	private boolean peakIsOpen() {
		// TODO Auto-generated method stub
		return false;
	}
	
}
