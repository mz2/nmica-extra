package net.derkholm.nmica.extra.app;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

import net.derkholm.nmica.apps.MotifScanner;
import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.model.MotifHitRecord;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;
import net.derkholm.nmica.model.motif.extra.BucketComparisonElement;
import net.derkholm.nmica.model.motif.extra.ScoredString;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;

@App(overview = "Determine the score cutoff for sequence motifs", generateStub = true)
@NMExtraApp(launchName = "nmcutoff")
public class MotifSetCutoffAssigner {

	/* required */
	private File motifFile;
	private File seqFile;
	private File backgroundModelFile;

	/* optional */
	private double minThreshold = -10.0;
	private File histogramOutputFile = null;
	private double bucketSize = 1.0;
	
	private File outputMotifFile = null; // when null, will output to original
	// file
	private double confidenceThreshold = 0.05;
	private HashSequenceDB sequences;
	private Hashtable<Motif, List<MotifHitRecord>> motifHitMap = new Hashtable<Motif, List<MotifHitRecord>>();
	private Hashtable<Motif, List<ScoredString>> enumSeqMap = new Hashtable<Motif,List<ScoredString>>();
	private MosaicSequenceBackground backgroundModel;
	
	@Option(help = "Input motifs")
	public void setMotifs(File motifFile) {
		this.motifFile = motifFile;
	}

	@Option(help = "Input sequences")
	public void setSeqs(Reader r) throws Exception {
		sequences = new HashSequenceDB();
		SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(r));
		while (si.hasNext()) {
			sequences.addSequence(si.nextSequence());
		}
	}

	@Option(help = "Suboptimal score threshold (default=-10.0)", userLevel = UserLevel.DEBUG, optional = true)
	public void setMinThreshold(double threshold) {
		this.minThreshold = threshold;
	}

	@Option(help="Background model in the NMICA background model XML format")
	public void setBackgroundModel(File f) throws Exception {
		XMLInputFactory factory = XMLInputFactory.newInstance();
        try {
    		XMLStreamReader r = factory.createXMLStreamReader(new FileReader(f));
    		Mosaic m = MosaicIO.readMosaic(r);
    		this.backgroundModel = new MosaicSequenceBackground(m.getDistributions(), m.getTransition());
    		r.close();
    	} catch (Exception ex) {
    		// ex.printStackTrace();
    		throw new Exception("Error loading background model. " +
    				"If you are using a background model created with " +
    				"an earlier version of NestedMICA, " +
    				"please try using the nmconvertbg program", ex);
    	}
	}

	@Option(help = "Output file for the score histogram", optional = true)
	public void setHistogramOutput(File histogramOutput) {
		this.histogramOutputFile = histogramOutput;
	}

	@Option(help = "Bucket size (default=true)", optional = true)
	public void setBucketSize(double d) {
		this.bucketSize = d;
	}

	@Option(help = "Confidence threshold (default=0.05)", optional = true)
	public void setConfThreshold(double d) {
		this.confidenceThreshold = d;
	}

	@Option(help = "Output file. If not specified, will override cutoff in the original file.", optional = true)
	public void setOut(File f) {
		this.outputMotifFile = f;
	}

	public void main(String[] args) throws Exception {
		Motif[] motifs = null;
		try {
			motifs = MotifIOTools
					.loadMotifSetXML(new FileInputStream(motifFile));
		} catch (FileNotFoundException e) {
			System.err.printf("Motif file not found: %s%n", motifFile
					.getAbsolutePath());
			System.exit(1);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.err.printf("Scanning motifs from against sequences...%n");
		MotifScanner scanner = new MotifScanner();
		scanner.setStoreHits(true);
		scanner.setScoreThreshold(minThreshold);
		scanner.scan(sequences, motifs);
		List<MotifHitRecord> hitRecords = scanner.hitRecords();
		for (Motif m : motifs) {
			motifHitMap.put(m, new ArrayList<MotifHitRecord>());
			for (MotifHitRecord rec : hitRecords) {
				if (rec.getMotif() == m) {
					motifHitMap.get(m).add(rec);					
				}
			}
		}
		
		hitRecords = null;
		scanner = null;
		
		System.err.printf(
			"Enumerating sequence words from background model " +
			"and weighting them according to the background model...%n");
		
		MotifMatchEnumerator motifEnumerator = new MotifMatchEnumerator();
		WordWeighter weighter = new WordWeighter();
		for (Motif m : motifs) {
			motifEnumerator.enumerateMatches(m, minThreshold);
			List<ScoredString> enums = motifEnumerator.storedHits();
			weighter.setBackgroundScoreForScoredStrings(enums, backgroundModel);
			enumSeqMap.put(m, enums);
		}
		motifEnumerator = null;
		weighter = null;
		
		System.err.println("Comparing expected and observed score histograms...");
		MotifMatchHistogramComparitor histogramComparitor = new MotifMatchHistogramComparitor();
		histogramComparitor.setBucketSize(bucketSize);
		histogramComparitor.setConfidence(confidenceThreshold);
		
		for (Motif m : motifs) {
			List<MotifHitRecord> realHits = motifHitMap.get(m);
			List<ScoredString> expHits = enumSeqMap.get(m);
			
			histogramComparitor.setConfidence(confidenceThreshold);
			histogramComparitor.setReal((List)realHits);
			histogramComparitor.setReference((List)expHits);
			
			double cutoff = histogramComparitor
								.determineCutoff(confidenceThreshold);
			System.err.printf("Cutoff for %s:%f%n",m.getName(),cutoff);
			m.setThreshold(cutoff);
		}
	}
}
