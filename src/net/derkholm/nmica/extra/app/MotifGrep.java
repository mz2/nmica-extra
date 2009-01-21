package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.motif.NMWeightMatrix;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;


@App(overview="A tool for outputting subsets of motif sets based on name lists or a regular expression", generateStub=true)
@NMExtraApp(launchName="nmgrep", vm=VirtualMachine.SERVER)
public class MotifGrep {

	private File motifs;
	private File list;
	private Pattern pattern;
	private String[] names;
	private String replaceWithStr = null;
	private boolean substring;
	private String prefix;
	private String annotationKey;
	private boolean matchSpecies;
	private boolean matchDesc;
	private boolean matchName;
	private int stripColumnsFromLeft;
	private int stripColumnsFromRight;
	
	@Option(help="List of motif names to match against. " +
			"Note that this is done by exact comparison, " +
			"not by finding matching substrings",optional=true)
	public void setList(File list) {
		this.list = list;
	}
	
	@Option(help="Strip specified number of columns from the left",optional=true)
	public void setStripColumnsFromLeft(int i) {
		stripColumnsFromLeft = i;
	}
	
	@Option(help="Strip specified number of columns from the right",optional=true)
	public void setStripColumnsFromRight(int i) {
		stripColumnsFromRight = i;
	}
	
	@Option(help="Allowed motif names (separated by spaces)",optional=true)
	public void setNames(String[] names) {
		this.names = names;
	}
	
	@Option(help="Find motifs whose name contains the specified string as a substring, " +
			"rather than looking for exact matches " +
			"(this switch works with -list and -names, default=false)", optional=true)
	public void setMatchSubstring(boolean b) {
		this.substring = b;
	}

	@Option(help="Regular expression to match against",optional=true)
	public void setExp(String str) {
		this.pattern = Pattern.compile(str);
	}
	
	@Option(help="Include the motif name in regular expression matching (default=true)",optional=true)
	public void setMatchName(boolean bool) {
		this.matchName = bool;
	}
	
	@Option(help="Include the species name in regular expression matching (default=false)",optional=true)
	public void setMatchSpecies(boolean bool) {
		this.matchSpecies = bool;
	}

	@Option(help="Include the description annotations in regular expression matching (default=false)",optional=true)
	public void setMatchDesc(boolean bool) {
		this.matchDesc = bool;
	}
	
	@Option(help="Print out the value for the specified annotation key",optional=true)
	public void setAnnotation(String str) {
		this.annotationKey = str;
	}
	
	@Option(help="Replacement string for the regular expression specified with -exp " +
			"(replaces all instances)",
			optional=true)
	public void setReplaceWith(String str) {
		this.replaceWithStr = str;
	}
	
	@Option(help="Add prefix to the motif names", optional=true) 
	public void setAddPrefix(String str) {
		this.prefix = str;
	}
	
	@Option(help="Input motifset file")
	public void setMotifs(File motifs) {
		this.motifs = motifs;
	}

	/**
	 * @param args
	 */
	public void main(String[] args) 
		throws Exception
	{
		Motif[] motifs = MotifIOTools.loadMotifSetXML(new FileReader(this.motifs));
		List<Motif> om = new ArrayList<Motif>();
		
		if (names != null && list != null) {
			System.err.println("Supply only either -names or -list");
			System.exit(1);
		}
		
		if (names != null) {
			for (String str : names) {
				for (Motif m : motifs) {
					if (substring ? m.getName().contains(str) : m.getName().equals(str)) {
						om.add(m);
						break;
					}
				}
			}
		} else if (list != null) {
			BufferedReader br = new BufferedReader(new FileReader(list));
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				for (Motif m : motifs) {
					if (substring ? m.getName().contains(line) : m.getName().equals(line)) {
						om.add(m);
						break;
					}
				}
			}
		} else if (pattern != null) {
			for (Motif m : motifs) {
				if (replaceWithStr != null) {
					Matcher matcher = pattern.matcher(m.getName());
					if (matcher.find()) {
						m.setName(matcher.replaceAll(replaceWithStr));
					}
					om.add(m);
				} else {
					if (matchName && (pattern.matcher(m.getName()).find())) 
					{om.add(m); continue;}
					if (matchDesc 
							&& m.getAnnotation().containsProperty("description") 
							&& (pattern.matcher(
								(CharSequence)m.getAnnotation()
									.getProperty("description")).find())) 
					{om.add(m); continue;}
					if (matchSpecies
							&& m.getAnnotation().containsProperty("species") 
							&& (pattern.matcher(
								(CharSequence)m.getAnnotation()
									.getProperty("species")).find())) 
					{om.add(m); continue;}
				}
			}
		} else {
			//just add every motif from input to the output set
			for (Motif m : motifs) om.add(m);
		}
		
		if (prefix != null) {
			for (Motif m : om) {
				m.setName(prefix + m.getName());
			}
		}
		
		if (annotationKey != null) {
			for (Motif m : om) {
				if (m.getAnnotation().containsProperty(annotationKey)) {
					System.out.printf("%s\t%s\t%s%n",
							m.getName(),
							annotationKey,
							m.getAnnotation().getProperty(annotationKey));
				}
			}
		} else {
			if (stripColumnsFromLeft > 0) {
				
				for (int m = 0; m < om.size(); m++) {
					Motif mot =  om.get(m);
					mot.setWeightMatrix(
						stripColumnsFromLeft(
								mot.getWeightMatrix(),
								stripColumnsFromLeft));
				}
			}
			
			if (stripColumnsFromRight > 0) {
				
				for (int m = 0; m < om.size(); m++) {
					Motif mot =  om.get(m);
					mot.setWeightMatrix(
						stripColumnsFromRight(
								mot.getWeightMatrix(),
								stripColumnsFromRight));
				}
			}
			
			MotifIOTools.writeMotifSetXML(
					System.out, 
					om.toArray(new Motif[0]));
		}
	}
	
	public WeightMatrix stripColumnsFromLeft(WeightMatrix inputWM, int count) throws IllegalAlphabetException {
		int colCount = inputWM.columns();
		
		if (colCount - count <= 0) {
			throw new IllegalArgumentException(
					"The input weight matrix hass too few columns.");
		}
		
		Distribution[] dists = new Distribution[colCount - count];
		
		for (int i = count-1; i < colCount; i++) {
			dists[i - count] = new SimpleDistribution(inputWM.getColumn(i));
		}
		
		return new NMWeightMatrix(dists,colCount - count, 0);
	}
	
	public WeightMatrix stripColumnsFromRight(WeightMatrix inputWM, int count) throws IllegalAlphabetException {
		int colCount = inputWM.columns();
		
		if (colCount - count <= 0) {
			throw new IllegalArgumentException(
					"The input weight matrix hass too few columns.");
		}
		
		Distribution[] dists = new Distribution[colCount - count];
		
		for (int i = 0; i < colCount - count; i++) {
			dists[i] = new SimpleDistribution(inputWM.getColumn(i+count));
		}
		
		return new NMWeightMatrix(dists,inputWM.columns() - count, 0);
	}
}