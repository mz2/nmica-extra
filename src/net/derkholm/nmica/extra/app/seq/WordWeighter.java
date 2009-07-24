package net.derkholm.nmica.extra.app.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;
import java.util.StringTokenizer;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.analysis.ScoredString;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ListTools;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Weight short sequences using a NestedMICA sequence background model", generateStub=true)
@NMExtraApp(launchName="nmweightwords", vm=VirtualMachine.SERVER)
public class WordWeighter {

	private MosaicSequenceBackground backgroundModel;
	private File wrdFile;

	@Option(help="Background model in the NMICA background model XML format")
	public void setBackgroundModel(File f) throws Exception {
		XMLInputFactory factory = XMLInputFactory.newInstance();
        try {
    		XMLStreamReader r = factory.createXMLStreamReader(new FileReader(f));
    		Mosaic m = MosaicIO.readMosaic(r);
    		backgroundModel = new MosaicSequenceBackground(m.getDistributions(), m.getTransition());
    		r.close();
    	} catch (Exception ex) {
    		// ex.printStackTrace();
    		throw new Exception("Error loading background model. " +
    				"If you are using a background model created with " +
    				"an earlier version of NestedMICA, " +
    				"please try using the nmconvertbg program", ex);
    	}
	}
	
	@Option(help="Sequence word file in the the output format of nmenumeratematches")
	public void setWords(File f) {
		this.wrdFile = f;
	}
	/**
	 * @param args
	 */
	public void main(String[] args) 
		throws Exception {
		Distribution[] mosaic = backgroundModel.getBackgroundDistributions();
		BufferedReader br = new BufferedReader(new FileReader(wrdFile));
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			StringTokenizer toke = new StringTokenizer(line);
			SymbolList word = DNATools.createDNA(toke.nextToken());
			double ww = 0;
			for (Distribution d : mosaic) {
				ww += weightWord(word, d);
			}
			System.out.println(line + "\t" + ww);
		}
	}
	
	public void setBackgroundScoreForScoredStrings(List<ScoredString> strings, MosaicSequenceBackground background) throws Exception {
		this.backgroundModel = background;
		Distribution[] mosaic = background.getBackgroundDistributions();
		
		for (ScoredString str : strings) {
			double bgscore = 0.0;
			for (Distribution d : mosaic) {
				bgscore += weightWord(DNATools.createDNA(str.getString()), d);
			}
			str.setBgScore(bgscore);
		}
	}

	public double weightWord(SymbolList word, Distribution dist)
		throws Exception {
		
		int order = dist.getAlphabet().getAlphabets().size();
		if (order != 2) {
			throw new Exception("Only 1st order models currently supported");
		} else {
			// Hard-coded order1 stuff.
			double p = 1.0;
			Symbol memory = DNATools.n();
			for (int i = 1; i <= word.length(); ++i) {
				p *= dist.getWeight(dist.getAlphabet().getSymbol(new ListTools.Doublet(memory, word.symbolAt(i))));
				memory = word.symbolAt(i);
			}
			return p;
		}
	}
}