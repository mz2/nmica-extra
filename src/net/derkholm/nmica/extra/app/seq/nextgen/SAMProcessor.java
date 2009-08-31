package net.derkholm.nmica.extra.app.seq.nextgen;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

import javax.naming.OperationNotSupportedException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.SequenceDB;
import org.bjv2.util.cli.Option;
import org.bjv2.util.cli.UserLevel;


public abstract class SAMProcessor {

	public enum IterationType {
		ONE_BY_ONE,
		MOVING_WINDOW,
		WITH_FREQUENCY,
		MAPPED_TO_REF
	}
	
	public enum QueryType {
		RECORD,
		CONTAINED,
		OVERLAP
	}
	
	protected SAMFileReader inReader;
	private SequenceDB seqDB = new HashSequenceDB();
	protected int qualityCutoff = 0;
	
	protected File indexFile;
	protected String in = "-";
	protected Map<String,Integer> refSeqLengths = new HashMap<String,Integer>();
	protected Map<String,Integer> readCounts = new HashMap<String,Integer>();

	protected int windowSize = 1;
	protected int frequency = 1;
	private IterationType iterationType = IterationType.ONE_BY_ONE;
	private QueryType queryType = QueryType.RECORD;
	private ArrayList<String> nameList;
	private boolean includeUnmapped = false;
	protected int extendedLength;
	protected int readLength;
	private boolean readLengthWasSet;
	private String currentRefSeqName;
	private int readQualityCutoff;


	@Option(help="Input reads (SAM/BAM formatted). Read from stdin if not specified.", optional=true)
	public void setMap(String in) {
		this.in = in;
	}
	
	public void setIterationType(IterationType type) {
		this.iterationType = type;
	}
	
	public void setQueryType(QueryType type) {
		this.queryType  = type;
	}
	
	@Option(help="Sequence window size around the current position (default=1)", optional=true)
	public void setWindowSize(int i) {
		this.windowSize = i;
	}
	
	@Option(help="Frequency of sequence windows (default=1)", optional=true)
	public void setWindowFreq(int i) {
		this.frequency = i;
	}
	
	@Option(help="Index file for the reads", optional=true)
	public void setIndex(File f) {
		this.indexFile = f;
	}
	
	@Option(help="Extend reads by specified number of nucleotides (bound by reference sequence ends)", optional=true)
	public void setExtendTo(int i) {
		this.extendedLength = i;
	}

	@Option(help="Reference sequence names and lengths in a TSV formatted file")
	public void setRefLengths(File f) throws NoSuchElementException, BioException, NumberFormatException, IOException {
		this.refSeqLengths = SAMProcessor.parseRefLengths(f);
		nameList = new ArrayList<String>(refSeqLengths.keySet());
		Collections.sort(
				nameList, 
				new Comparator<String>() {
					public int compare(String str1, String str2) {
						return str1.compareTo(str2);
					}
				});
	}
	
	public static Map<String, Integer> parseRefLengths(File f) {
		Map<String, Integer> refSeqLengths = new HashMap<String, Integer>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(f));
			
			String line = null;
			while ((line = reader.readLine()) != null) {
				StringTokenizer tok = new StringTokenizer(line,"\t");
				refSeqLengths.put(tok.nextToken(), Integer.parseInt(tok.nextToken()));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return refSeqLengths;
	}

	@Option(
		help="Read counts in a TSV formatted file " +
				"(reference seq name + tab + counts for that ref seq)",
				optional=true)
	public void setReadCounts(File f) {
		this.readCounts = SAMProcessor.parseReadCounts(f);
	}
	
	static Map<String, Integer> parseReadCounts(File f) {
		Map<String,Integer> readCounts = new HashMap<String,Integer>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(f));
			
			String line = null;
			while ((line = reader.readLine()) != null) {
				StringTokenizer tok = new StringTokenizer(line, "\t");
				readCounts.put(tok.nextToken(), Integer.parseInt(tok.nextToken()));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return readCounts;	
	}

	@Option(help="Sequencing read length", optional=true)
	public void setReadLength(int i) {
		this.readLengthWasSet = true;
		this.readLength = i;
	}
	
	public void setMappingQualityAbove(int quality) {
		this.qualityCutoff = quality;
	}
	
	public void setReadQualityAbove(int quality) {
		this.readQualityCutoff = quality;
	}
	
	@Option(help="Include unmapped reads (default=false)", optional=true, userLevel=UserLevel.DEBUG)
	public void setIncludeUnmapped(boolean b) {
		this.includeUnmapped  = b;
	}
	
	//override in subclass to handle QueryType.CONTAINED and QueryType.OVERLAP
	public void process(
			List<SAMRecord> recs, 
			String refName, 
			int begin, 
			int end, 
			int seqLength) {
	
	}
	
	//override in subclass to handle QueryType.ONE_BY_ONE
	public void process(SAMRecord rec, int readIndex) {
		
	}
	
	//override in subclass
	public void processAndClose(CloseableIterator<SAMRecord> recs, String refName, int len) {
	
	}
	
	public void initializeSAMReader() {
		if (this.in.equals("-")) {
			if (queryType != QueryType.RECORD) {
				System.err.println("Query type -record is the only allowed query type when reading from stdin");
				System.exit(1);
			}
			this.inReader = new SAMFileReader(System.in);
		} else {
			if (indexFile == null && (this.queryType.equals(QueryType.CONTAINED) || this.queryType.equals(QueryType.OVERLAP))) {
				System.err.println("Index file was not specified but is required for query types 'contained' and 'overlap'");
				System.exit(2);
			}
			
			if (indexFile != null) {
				this.inReader = new SAMFileReader(new File(in),indexFile);				
			} else {
				this.inReader = new SAMFileReader(new File(in));
			}
		}
		this.inReader.setValidationStringency(ValidationStringency.SILENT);
	}

	//the main method in subclasses can pretty much look like this (you can customise of course)
	public void main(String[] args) throws Exception {
		initializeSAMReader();
		process();
	}
	
	public void process() throws BioException {
		int halfWindow = Math.max(1,this.windowSize / 2);

		if (iterationType == IterationType.ONE_BY_ONE) {
			int excludedReads = 0;
			int readCount = 0;

			for (SAMRecord record : inReader) {
				if ((readCount++ % frequency) != 0) continue;
				
				int quality = record.getMappingQuality();
				if ((qualityCutoff > 0) && (quality < qualityCutoff)) {
					excludedReads += 1;
					continue;
				}
				
				if (record.getReadUnmappedFlag()) {
					excludedReads += 1;
					continue;
				}
				
				if ((readQualityCutoff > 0) && (record.getIntegerAttribute("UQ") < readQualityCutoff)) {
					excludedReads += 1;
					continue;
				}
				
				if (this.extendedLength > 0) {
					ExtendReads.extendReadBy(
							record, 
							refSeqLengths, 
							this.extendedLength - this.readLength);
				}
				process(record, readCount);
			}
			System.err.printf(
				"Excluded %d reads (%.2f%%)%n", 
				excludedReads, 
				(double)excludedReads / (double)readCount * 100.0);			
		} else if (iterationType == IterationType.MOVING_WINDOW) {
			
			final List<SAMRecord> recs = new ArrayList<SAMRecord>();
			
			for (String seqName : nameList) {
				System.err.printf("Processing %s%n",seqName);
				
				int windowCenter = halfWindow;
				int len = refSeqLengths.get(seqName);
				while ((windowCenter + halfWindow) < len) {
					CloseableIterator<SAMRecord> recIterator;
					
					if (this.extendedLength > 0) {
						int extendedStart = windowCenter - extendedLength;
						int extendedEnd = windowCenter + extendedLength;
						
						recIterator = this.query(seqName, extendedStart, extendedEnd);
						iterateAndFilterToList(recIterator,windowCenter,recs);
					} else {
						recIterator = this.query(seqName, windowCenter - halfWindow, windowCenter + Math.max(1,halfWindow));
						iterateAndFilterToList(recIterator,windowCenter,recs);
					}
					recIterator.close();
					
					process(recs,seqName,windowCenter - halfWindow,windowCenter + Math.max(1,halfWindow),len);
					recs.clear();
					windowCenter += frequency;
				}
				System.err.printf(".");
			}
			
		} else if (iterationType == IterationType.WITH_FREQUENCY) {
			windowSize = frequency;
			
			final List<SAMRecord> recs = new ArrayList<SAMRecord>();
			for (String seqName : nameList) {
				System.err.printf("Processing %s%n",seqName);
				
				int windowBegin = 0;
				int len = refSeqLengths.get(seqName);
				
				while ((windowBegin + frequency)  < len) {
					CloseableIterator<SAMRecord> recIterator = this.query(seqName, windowBegin, windowBegin + frequency);
					iterateAndFilterToList(recIterator, windowBegin + (frequency / 2), recs);
					recIterator.close();
					
					process(recs,seqName,windowBegin,windowBegin + frequency,len);
					recs.clear();
				}
				windowBegin += frequency;
			}
		} else if (iterationType == IterationType.MAPPED_TO_REF) {
			for (String seqName : nameList) {
				System.err.printf("Processing %s%n",seqName);
				
				this.setCurrentRefSeqName(seqName);
				CloseableIterator<SAMRecord> recs = 
					inReader.queryContained(seqName, 0, refSeqLengths.get(seqName));
				processAndClose(recs, seqName, refSeqLengths.get(seqName));
			}
		}
	}

	protected void setCurrentRefSeqName(String seqName) {
		this.currentRefSeqName = seqName;
	}

	private void iterateAndFilterToList(
			CloseableIterator<SAMRecord> recIterator,
			int windowCenter,
			final List<SAMRecord> recs) {
		int winStart = windowCenter - (windowSize / 2);
		int winEnd = windowCenter + (windowSize / 2);
		
		while (recIterator.hasNext()) {
			SAMRecord rec = recIterator.next();
			if (rec.getMappingQuality() < this.qualityCutoff) continue;
			
			
			if (extendedLength > 0) {
				boolean onPositiveStrand = !rec.getReadNegativeStrandFlag();
				
				if (onPositiveStrand) {
					if ((rec.getAlignmentStart() + extendedLength) < winStart) {
						
						System.err.printf("%d - %d on + strand cannot be extended to hit win start at %d%n", 
								rec.getAlignmentStart(), 
								rec.getAlignmentStart() + extendedLength,
								winStart);
						continue; // if you can't extend the read to hit the win start pos
					}
					if (rec.getAlignmentStart() > winEnd) {

						System.err.printf("%d - %d on + strand cannot be extended to hit win end at %d%n", 
								rec.getAlignmentStart(), 
								rec.getAlignmentStart() + extendedLength,
								winEnd);
						continue; // if the read doesn't start before the window ends
					}
				} else {
					if (rec.getAlignmentEnd() < winStart) {
						System.err.printf("%d - %d on - strand cannot be extended to hit win start at %d%n", 
								rec.getAlignmentStart(), 
								rec.getAlignmentStart() + extendedLength);
						continue; //if the read doesn't start before the window starts
					}
					if ((rec.getAlignmentEnd() - extendedLength) > winEnd) {
						System.err.printf("%d - %d on - strand cannot be extended to hit win end at %d%n", 
								rec.getAlignmentStart(), 
								rec.getAlignmentStart() + extendedLength,
								winEnd);
						continue; // if the read can't be extended to hit inside the window (min coordinate smaller than window's end)
					}
				}
				System.err.printf("YES! Read %d - %d (%s strand) overlaps with the target range %d - %d%n",
						rec.getAlignmentStart(), 
						rec.getAlignmentStart() + extendedLength, 
						rec.getReadNegativeStrandFlag() ? "-" : "+",
						winStart, winEnd);
			}
			
			recs.add(rec);
		}
		return;
	}
	
	private CloseableIterator<SAMRecord> query(String seqName, int begin, int end) {
		if (this.queryType == QueryType.CONTAINED) {
			return inReader.queryContained(seqName, begin, end);			
		} else {
			return inReader.queryOverlapping(seqName, begin, end);
		}
	}

	public void setOut(File f) {
		// TODO Auto-generated method stub
		return;
	}
}
