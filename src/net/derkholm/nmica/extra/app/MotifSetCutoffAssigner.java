package net.derkholm.nmica.extra.app;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

import net.derkholm.nmica.apps.MotifScanner;
import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.model.MotifHitRecord;
import net.derkholm.nmica.model.analysis.ScoredString;
import net.derkholm.nmica.model.motif.Mosaic;
import net.derkholm.nmica.model.motif.MosaicIO;
import net.derkholm.nmica.model.motif.MosaicSequenceBackground;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifIOTools;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.utils.ChangeVetoException;
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
	private File sequences;
	private ConcurrentHashMap<Motif, List<MotifHitRecord>> motifHitMap = new ConcurrentHashMap<Motif, List<MotifHitRecord>>();
	private ConcurrentHashMap<Motif, List<ScoredString>> enumSeqMap = new ConcurrentHashMap<Motif,List<ScoredString>>();
	private MosaicSequenceBackground backgroundModel;
	private double defaultThreshold;
	private int threads = 1;
	private ExecutorService threadPool;
	
	@Option(help = "Input motifs")
	public void setMotifs(File motifFile) {
		this.motifFile = motifFile;
	}

	@Option(help = "Input sequences")
	public void setSeqs(File f) throws Exception {
		sequences = f;
	}

	@Option(help = "Suboptimal score threshold (default=-10.0)", userLevel = UserLevel.DEBUG, optional = true)
	public void setSuboptimalThreshold(double threshold) {
		this.minThreshold = threshold;
	}
	
	@Option(help = "Default score threshold applied to motifs for which threshold could not be determined " +
			"(an annotation 'default_threshold_used' will also be added to these cases)", optional = true)
	public void setDefaultThreshold(double d) {
		this.defaultThreshold = d;
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
	
	@Option(help="Number of threads (default = 1)", optional = true)
	public void setThreads(int threads) {
		if (threads < 1) {
			System.err.println("-threads needs to be >= 1");
			System.exit(1);
		}
		this.threads  = threads;
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

		if (this.threads > motifs.length) {
			System.err.println("Specified number of threads is larger than the number of motifs in the input. " +
					"Will use only " + motifs.length + " threads.");
		}
		this.threadPool = Executors.newFixedThreadPool(this.threads);
		

		System.err.printf("Scanning motifs against sequences...%n");
		
		List<Future<Boolean>> scanFutures = new ArrayList<Future<Boolean>>();
		for (Motif m : motifs) {
			
			scanFutures.add(threadPool.submit(
					new ScanTask(
							m, 
							sequences, 
							motifHitMap,minThreshold)));
		}
		for (Future<Boolean> sf : scanFutures) {
			if (!sf.get().booleanValue()) {
				System.err.println("Unexpected failure: one of the scanning tasks failed.");
				System.exit(2);
			}
		}
		threadPool.shutdown();
		
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
		this.threadPool = Executors.newFixedThreadPool(this.threads);
		List<Future<Motif>> motifFutures = new ArrayList<Future<Motif>>();
		for (Motif m : motifs) {
			threadPool.submit(new CutoffTask(
					m, 
					motifHitMap.get(m), 
					enumSeqMap.get(m), 
					bucketSize, 
					confidenceThreshold, 
					defaultThreshold));
		}
		
		for (Future<Motif> mf : motifFutures) {
			Motif m = mf.get();
			System.err.printf("%s\t%f%n",m.getName(),m.getThreshold());
		}
		threadPool.shutdown();
		
		OutputStream os = null;
		if (outputMotifFile == null) {
			os = System.out;
		} else {
			os = new BufferedOutputStream(new FileOutputStream(outputMotifFile));
		}
		
		MotifIOTools.writeMotifSetXML(os, motifs);
	}
	
	private static class ScanTask implements Callable<Boolean> {
		private final ConcurrentHashMap<Motif,List<MotifHitRecord>>  motifHitMap;
		private File seqFile;
		private HashSequenceDB sequences;
		private final Motif motif;
		private double minThreshold;
		
		public ScanTask(
				Motif motif, 
				File seqFile,
				ConcurrentHashMap<Motif,List<MotifHitRecord>> motifHitMap,
				double minThreshold) throws FileNotFoundException, ChangeVetoException, NoSuchElementException, BioException {
			this.motif = motif;
			this.seqFile = seqFile;
			
			HashSequenceDB seqDB = sequences = new HashSequenceDB();
			SequenceIterator si = SeqIOTools.readFastaDNA(new BufferedReader(new FileReader(seqFile)));
			while (si.hasNext()) {sequences.addSequence(si.nextSequence());}
			
			this.motifHitMap = motifHitMap;
			this.minThreshold = minThreshold;
		}
		
		public Boolean call() {
			MotifScanner scanner = new MotifScanner();
			scanner.setStoreHits(true);
			scanner.setScoreThreshold(minThreshold);
			scanner.setFormat(MotifScanner.Format.GFF);
			try {
				System.err.printf("Scanning sequences against %s...%n", motif.getName());
				scanner.scan(sequences, new Motif[]{motif});
			} catch (Exception e) {
				throw new BioError("Scanning failed for motif " + this.motif);
			}
			List<MotifHitRecord> hitRecords = scanner.hitRecords();
			
			for (MotifHitRecord rec : hitRecords) {
				System.err.println(rec.getScore());
				if (rec.getMotif() == motif) {
					hitRecords.add(rec);					
				}
			}
			motifHitMap.put(motif, hitRecords);
			
			return new Boolean(true);
		}
	}
	
	private static class CutoffTask implements Callable<Motif> {
		private final MotifMatchHistogramComparitor histogramComparitor;
		private final Motif motif;
		private final List<MotifHitRecord> realHits;
		private final List<ScoredString> expHits;
		private final double confidenceThreshold;
		private final double defaultThreshold;
		
		public CutoffTask(
				Motif m,
				List<MotifHitRecord> realHits,
				List<ScoredString> expHits,
				double bucketSize, 
				double confidenceThreshold, 
				double defaultThreshold) {
			this.motif = m;
			this.realHits = realHits;
			this.expHits = expHits;
			this.confidenceThreshold = confidenceThreshold;
			this.defaultThreshold = defaultThreshold;
			
			this.histogramComparitor = new MotifMatchHistogramComparitor();
			histogramComparitor.setBucketSize(bucketSize);
			histogramComparitor.setConfidence(confidenceThreshold);
			
			histogramComparitor.setConfidence(confidenceThreshold);
			histogramComparitor.setReal((List)realHits);
			histogramComparitor.setReference((List)expHits);
		}

		public Motif call() throws Exception {
			double cutoff = histogramComparitor.determineCutoff(confidenceThreshold);

			if (Double.isNaN(cutoff)) {
				System.err.printf("Setting default cutoff %f for motif %s%n",cutoff,motif.getName());
				motif.setThreshold(this.defaultThreshold);
				motif.getAnnotation().setProperty("default_threshold_used", "" + this.defaultThreshold);
			}
			motif.setThreshold(cutoff);
			return motif;
		}
	}
	
}
