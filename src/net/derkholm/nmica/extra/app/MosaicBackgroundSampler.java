package net.derkholm.nmica.extra.app;

import java.io.InputStream;
import java.util.Iterator;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;

import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Output sequences sampled from a mosaic background", generateStub=true)
@NMExtraApp(launchName="nmbgsample", vm=VirtualMachine.SERVER)
public class MosaicBackgroundSampler {
	
	private MosaicSequenceBackground backgroundModel;
	
	private int count = 1;
	private int length = 100;

	@Option(help="Sampled sequence count (default=1)", optional=true)
	public void setCount(int i) {
		this.count = i;
	}
	
	@Option(help="Sampled sequence length (default=100)", optional=true) 
	public void setLength(int len) {
		this.length = len;
	}
	
	@Option(help="The background model to use", optional=true)
    public void setBackgroundModel(InputStream is)
        throws Exception {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        
        try {
    		XMLStreamReader r = factory.createXMLStreamReader(is);
    		Mosaic m = MosaicIO.readMosaic(r);
    		backgroundModel = new MosaicSequenceBackground(m.getDistributions(), m.getTransition());
    		r.close();
    	} catch (Exception ex) {
    		// ex.printStackTrace();
    		throw 
    			new Exception(
					"Error loading background model." +
					" If you are using a background model created with an earlier version of NestedMICA, " +
					"please try using the nmconvertbg program", ex);
    	}
    }
	
	public void main(String[] args) throws Exception {
		DP dp = backgroundModel.getBackgroundDP();

		if ((length % backgroundModel.getMosaicOrder()) != 0) {
			System.out.println("Sequence length has to be divisible " +
					"by the background model order (" + backgroundModel.getMosaicOrder() + ")");
			System.exit(1);
		}
		length = length / backgroundModel.getMosaicOrder();
		
		for (int c = 0; c < count; c++) {
			StatePath path = dp.generate(length);
			SymbolList symList = path.symbolListForLabel(StatePath.SEQUENCE);

			List states = path.toList();
			
			Alphabet alphabet = null;
			
			if (path.getAlphabet().getName().contains("DNA"))
				alphabet = DNATools.getDNA();
			else if (path.getAlphabet().getName().contains("protein"))
				alphabet = ProteinTools.getAlphabet();
			else {
				System.out.println("Unsupported alphabet: " + alphabet.getName());
				System.exit(1);
			}
			
			SymbolTokenization symTok = alphabet.getTokenization("token");
			
			System.out.println(">seq" + c);
			for (Iterator it = symList.iterator(); it.hasNext();) {
				Object o = it.next();
				AtomicSymbol sym = (AtomicSymbol)o;
				List syms = sym.getSymbols();
				
				for (int i = 0; i < backgroundModel.getMosaicOrder(); i++) {
					Object symO = syms.get(i);
					Symbol bsym = (Symbol) symO;
					String symN = symTok.tokenizeSymbol(bsym);
					System.out.print(symN);
				}
			}
			System.out.println();
		}
	}
}
