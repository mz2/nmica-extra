package net.derkholm.nmica.extra.app;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview = "Get noncoding sequences from Ensembl for motif discovery", 
		generateStub = true)
@NMExtraApp(launchName = "nmgetensemblseq", 
			vm = VirtualMachine.SERVER)
public class RetrieveEnsemblSequences {
	
	/* Sequence region processing */
	private boolean repeatMask = true;
	private boolean excludeTranslations = true;
	private int featherTranslationsBy = 0;
	private int featherRegionsBy = 0;

	/* Sequence filtering properties 
	 * id = list of identifiers
	 * idType = type of identifier ('ensembl_gene' and 'stable_id' for now)*/
	private String[] ids;
	private String idType = "ensembl_gene";

	/* Sequence region properties */
	private int threePrimeBegin,threePrimeEnd;
	private int fivePrimeBegin,fivePrimeEnd;
	
	/* Ensembl settings */
	private String username = "anonymous";
	private String password = "";
	private String host = "ensembldb.ensembl.org";
	private String database;
	private int port = 5306;
	private int schemaVersion = 54;
	
	@Option(help="Repeat mask the sequences (default=true)", optional=true)
	public void setRepeatMask(boolean b) {
		this.repeatMask = b;
	}
	
	@Option(help="Exclude translations (default=true)", optional=true)
	public void setExcludeTranslations(boolean b) {
		this.excludeTranslations  = b;
	}
	
	@Option(help="Feather translations by the specified number of nucleotides", optional=true)
	public void setFeatherTranslationsBy(int i) {
		this.featherTranslationsBy = i;
	}
	
	@Option(help="Feather retrieved sequence regions by the specified number of nucleotides", optional = true)
	public void setFeatherRegionsBy(int i) {
		this.featherRegionsBy = i;
	}
	
	@Option(help="Filter returned sequence set based on identifiers")
	public void setFilterById(String[] ids) {
		this.ids = ids;
	}
	
	@Option(help="Identifier type to use when filtering sequences identifiers (default = ensembl_gene, possible values: ensembl_gene|stable_id)")
	public void setIdType(String idType) {
		this.idType = idType;
	}
	
	@Option(help="Get three prime UTR sequences. " +
			"Example: '-threePrimeUTR -200 200' gets you sequence regions " +
			"from 200bp upstream to 200bp downstream of transcription " +
			"start sites.", optional=true)
	public void setThreePrimeUTR(String[] coordStrs) {
		if (coordStrs.length != 2) {
			System.err.printf("-threePrimeUTR requires two arguments: begin and end coordinate (%d given)",coordStrs.length);
			System.exit(1);
		}
		
		this.threePrimeBegin = Integer.parseInt(coordStrs[0]);
		this.threePrimeEnd = Integer.parseInt(coordStrs[1]);
	}
	
	@Option(help="Get five prime UTR sequences. " +
			"Example: '-threePrimeUTR -200 200' gets you sequence regions " +
			"from 200bp upstream to 200bp downstream of transcription " +
			"start sites.", optional=true)
	public void setFivePrimeUTR(String[] coordStrs) {
		if (coordStrs.length != 2) {
			System.err.printf("-fivePrimeUTR requires two arguments: begin and end coordinate (%d given)",coordStrs.length);
			System.exit(2);
		}
		
		this.fivePrimeBegin = Integer.parseInt(coordStrs[0]);
		this.fivePrimeEnd = Integer.parseInt(coordStrs[1]);
	}
	
	@Option(help="Ensembl username (default=anonymous)", optional=true)
	public void setUser(String user) {
		this.username = user;
	}
	
	@Option(help="Ensembl password (default='')", optional=true)
	public void setPassword(String password) {
		this.password = password;
	}
	
	@Option(help="Ensembl hostname (default=ensembldb.ensembl.org", optional=true)
	public void setHost(String hostname) {
		this.host = hostname;
	}
	
	@Option(help="Ensembl database port (default=5306)", optional=true)
	public void setPort(int port) {
		this.port = port;
	}
		
	@Option(help="Ensembl database to retrieve sequences from (e.g. 'danio_rerio_core_54_8')", optional=true)
	public void setDatabase(String db) {
		this.database = db;
	}
	
	@Option(help="Ensembl schema version (default = 54)", optional=true)
	public void setSchemaVersion(int ver) {
		this.schemaVersion = ver;
	}
	
	public void main(String[] argv) throws Exception {
		String dbURL = String.format("jdbc:mysql:/%s:%d/%s", this.host, this.port, this.database);
		//EnsemblConnection conn = new EnsemblConnection(dbURL, username, password, schemaVersion);
		
		System.out.println(">foo\nactacatcatgggacata\n");
		
		//sequence retrieval Ensembl magic goes here.
	}
}
