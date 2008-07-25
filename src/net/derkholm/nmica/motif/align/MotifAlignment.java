package net.derkholm.nmica.motif.align;

import java.beans.PropertyChangeSupport;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import net.derkholm.nmica.model.metamotif.Dirichlet;
import net.derkholm.nmica.model.metamotif.DirichletParamEstimator;
import net.derkholm.nmica.model.metamotif.MetaMotif;
import net.derkholm.nmica.motif.Motif;
import net.derkholm.nmica.motif.MotifComparisonMatrixBundle;
import net.derkholm.nmica.motif.MotifComparitorIFace;
import net.derkholm.nmica.seq.WmTools;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;


public class MotifAlignment implements 
	Iterable<MotifAlignment.MotifAlignmentElement> {
	protected List<MotifAlignmentElement> motifElems;
	protected SortedSet<MotifAlignmentElement> offsetSortedMotifs;
	protected HashMap<Motif, MotifAlignmentElement> motifsToAlignmentElems;
	protected SortedSet<Integer> alignedMotifIndices;
	protected Set<Motif> flippedMotifs;
	
	private final PropertyChangeSupport pcs = new PropertyChangeSupport(this);
	
	protected String name;
	protected final FiniteAlphabet alphabet;
	protected Motif[] allMotifs;
	
	
	private MotifAlignment(FiniteAlphabet alphabet) {
		this.alphabet = alphabet;
		
		motifElems = new ArrayList<MotifAlignmentElement>();
		offsetSortedMotifs = 
			new TreeSet<MotifAlignmentElement>(
					new MotifAlignmentElementOffsetComparator());
		
		motifsToAlignmentElems = 
			new HashMap<Motif, MotifAlignmentElement>();
		
		alignedMotifIndices = new TreeSet<Integer>();
		flippedMotifs = new HashSet<Motif>();
	}
	
	public MotifAlignment(Motif[] motifs, MotifComparitorIFace mc) 
	throws Exception {
		this((FiniteAlphabet) motifs[0].getWeightMatrix().getAlphabet());
		allMotifs = motifs; /* used for indexing */
		List<Motif> motifList = new ArrayList<Motif>(Arrays.asList(motifs));

		MotifComparisonMatrixBundle mb = mc.getComparisonMatrixWithOffsets(
				motifList.toArray(new Motif[motifList.size()]));
		
		MotifPairWithOffset mp = highestScoringPair(mb);
		mp.getM1();
		mp.getM2();
		
		addMotif(mp.getM1(), 0, false);
		addMotif(mp.getM2(), mp.getOffset(), mp.isFlipped());
		
		motifList.remove(mp.getM1());
		motifList.remove(mp.getM2());
		
		//System.err.printf(" # M1 : %s M2: %s offset: %d score: %.3g flipped: %b %n",
		//		mp.getM1().getName(), 
		//		mp.getM2().getName(),
		//		mp.getOffset(),
		//		mp.getScore(),
		//		mp.isFlipped());
		
		if (mp.isFlipped()) {
			String oldConsensus = this.consensus(mp.getM2());
			mp.getM2()
				.setWeightMatrix(
					WmTools.reverseComplement(
						mp.getM2().getWeightMatrix()));
			String newConsensus = this.consensus(mp.getM2());
			flippedMotifs.add(mp.getM2());
			//System.err.printf("Flipping %s (%s --> %s)%n", 
			//		mp.getM2().getName(), 
			//		oldConsensus, newConsensus);
		}
		
		//System.err.printf("%n%s",this.alignmentConsensusString());
		//System.err.print("Flipped: ");
		//for (Motif m : flippedMotifs) System.err.printf("%s ",m.getName());
		//System.err.println("\n");
		
		/* then iterate through the remaining ones */
		while (motifList.size() > 0) {
			
			mp = highestScoringPairWithAlreadyAligned(mb);
			
			
			if (mp.isFlipped() &! (flippedMotifs.contains(mp.getM1()))) {
				//if (flippedMotifs.contains(mp.getM1()) &&) 
				
				//if (flippedMotifs.contains(mp.getM1()) &)
				String oldConsensus = this.consensus(mp.getM2());
				flippedMotifs.add(mp.getM2());
				mp.getM2()
					.setWeightMatrix(
						WmTools.reverseComplement(
							mp.getM2().getWeightMatrix()));
				String newConsensus = this.consensus(mp.getM2());
				flippedMotifs.add(mp.getM2());
				System.err.printf(" > Flipping %s (%s --> %s)%n", 
						mp.getM2().getName(), 
						oldConsensus, newConsensus);
			} else if (!mp.isFlipped() && (flippedMotifs.contains(mp.getM1()))) {
				String oldConsensus = this.consensus(mp.getM2());
				flippedMotifs.add(mp.getM2());
				mp.getM2()
					.setWeightMatrix(
						WmTools.reverseComplement(
							mp.getM2().getWeightMatrix()));
				String newConsensus = this.consensus(mp.getM2());
				flippedMotifs.add(mp.getM2());
				System.err.printf(" >> Flipping %s (%s --> %s)%n", 
						mp.getM2().getName(), 
						oldConsensus, newConsensus);
			}
			
			//System.err.printf(" > M1 : %s M2: %s offset: %d score: %.3g flipped: %b %n",
			//		mp.getM1().getName(), 
			//		mp.getM2().getName(),
			//		mp.getOffset(),
			//		mp.getScore(),
			//		mp.isFlipped());
			
			
			int offset;
			if (mp.isFlipped()) {
				offset = offset(mp.getM1()) + mp.getOffset();
			} else {
				offset = offset(mp.getM1()) + mp.getOffset();
			}
			
			addMotif(mp.getM2(),
						offset,
						mp.isFlipped());
			
			//System.err.printf("%n%s",this.alignmentConsensusString());
			//System.err.print("Flipped: ");
			//for (Motif m : flippedMotifs) System.err.printf("%s ",m.getName());
			//System.err.println("\n");
			motifList.remove(mp.getM2());
		}
		
		//System.err.println();
		//System.err.print("Flipped: ");
		//for (Motif m : flippedMotifs) System.err.printf("%s ",m.getName());
		//System.err.println("\n");
	}

	public MotifAlignment(MotifAlignment motifAlignment) {
		this(motifAlignment.alphabet);
		this.setName(motifAlignment.getName());
		
		for (MotifAlignmentElement elem : motifAlignment.motifElems)
			addMotif(new Motif(elem.getMotif()), 
					 elem.getOffset(), 
					 elem.isFlipped());
		
	}

	public MotifAlignment alignmentWithZeroOffset() {
		int minOffset = minOffset();
		
		MotifAlignment alignment = new MotifAlignment(alphabet);
		alignment.setName(name);
		for (MotifAlignmentElement elem : this.motifElems) {
			alignment.addMotif(elem.getMotif(), elem.getOffset() - minOffset, elem.isFlipped());
		}
		
		return alignment;
	}
		
	public void addMotif(Motif m, int offset, boolean flipped) {
		//System.err.println("Adding " + m.getName());
		//System.err.println("+"+m.getName()+" offset:" + offset);
		if (!contains(m)) {
			MotifAlignmentElement me = new MotifAlignmentElement(m,offset,flipped);
			motifElems.add(me);
			motifsToAlignmentElems.put(m, me);
			offsetSortedMotifs.add(me);
			
			/* This aligned motif indexing is used 
			 * with the MotifAlignment(Motif[] motifs,MotifComparitorIFace mc) 
			 * constructor */
			if (allMotifs != null) {
				int i = -1;
				//System.out.println("All motifs " + allMotifs);
				for (int index = 0; i < allMotifs.length; index++)
					if (allMotifs[index] == m) { i = index; break; }
				
				if (i == -1)
					throw new IllegalStateException(
						"Motif " + m + 
						" is not contained in the motif array given");
				
				alignedMotifIndices.add(i);
			}
		} else {
			throw new IllegalArgumentException(
					"Motif " + m.getName() + 
					" is already present in the alignment (offset " + offset + ")");
		}
	}
	
	public boolean removeMotif(Motif m) {
		if (motifsToAlignmentElems.containsKey(m)) {
			MotifAlignmentElement me = motifsToAlignmentElems.get(m);
			motifElems.remove(me);
			offsetSortedMotifs.remove(me);
			motifsToAlignmentElems.remove(m);
			return true;
		} else {
			return false;
		}
	}
	
	public boolean contains(Motif m) {
		return motifElems.contains(motifsToAlignmentElems.get(m));
	}
	
	public int offset(Motif m) {
		//System.out.println("offset for " + m.getName());
		//System.out.println(motifsToAlignmentElems.get(m));
		return motifsToAlignmentElems.get(m).getOffset();
	}

	public int numberOfMotifs() {
		return motifElems.size();
	}
	
	public int totalLength() {
		int maxLength = 0;
		int minOffset = minOffset();
		
		for (MotifAlignmentElement elem : motifElems) {
			WeightMatrix wm = elem.getMotif().getWeightMatrix();
			int columns = wm.columns();
			int offset = elem.getOffset();
			int thisMaxLength = columns + offset - minOffset;
			if (thisMaxLength > maxLength)
				maxLength = thisMaxLength;
		}
		
		return maxLength;
	}

	public MotifAlignmentElement[] toArray(MotifAlignmentElement[] motifs) {
		return this.motifElems.toArray(motifs);
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	private int maxOffset() {
		return offsetSortedMotifs.last().getOffset();
	}

	private int minOffset() {
		return offsetSortedMotifs.first().getOffset();
	}
	

	public class MotifAlignmentElement {
		private Motif motif;
		private int offset;
		private boolean flipped;
		
		public Motif getMotif() {
			return motif;
		}

		public void setMotif(Motif motif) {
			this.motif = motif;
		}

		public int getOffset() {
			return offset;
		}

		public void setOffset(int offset) {
			this.offset = offset;
		}

		public MotifAlignmentElement(Motif m, int offset, boolean flipped) {
			this.motif = m;
			this.offset = offset;
			this.flipped = flipped;
		}
		
		public Distribution getPosition(int position) {
			int offsetPos = position - offset;
			
			if (motif.getWeightMatrix().columns() > offsetPos && offsetPos >= 0) {
				return motif.getWeightMatrix().getColumn(offsetPos);
			} else
				return null;
		}

		public boolean isFlipped() {
			return flipped;
		}
	}
	
	public class MotifAlignmentElementOffsetComparator 
		implements Comparator<MotifAlignmentElement> {

		public int compare(MotifAlignmentElement o1, MotifAlignmentElement o2) {
			if (o1.offset < o2.offset)
				return -1;
			else if (o1.offset > o2.offset)
				return 1;
			else return 0;
		}
	}
	
	public class MotifAlignmentElementAlphabeticComparator
		implements Comparator<MotifAlignmentElement> {

		public int compare(MotifAlignmentElement o1, MotifAlignmentElement o2) {
			return o1.motif.getName().compareTo(o2.motif.getName());
		}
		
	}

	private Motif[] _motifs() {
		Motif[] ms = new Motif[motifElems.size()];
		for (int i = 0; i < motifElems.size(); i++) {
			ms[i] = motifElems.get(i).getMotif();
			//ms[i].setName(ms[i].getName() + "_" + motifElems.get(i).getOffset());
		}
		
		return ms;
	}
	
	public Motif[] motifs() {
		writeOffsetAnnotations();
		return _motifs();
	}

	
	public Distribution[] distributionsAtPosition(int position) {
		List<Distribution> dists = new ArrayList<Distribution>();
		int i = 0;
		for (MotifAlignmentElement elem : motifElems) {
			Distribution d = elem.getPosition(position);
			if (d != null) dists.add(d);
		}
		
		return dists.toArray(new Distribution[dists.size()]);
	}
	
	public MetaMotif metamotif(boolean trimEnds) 
		throws IllegalSymbolException, IllegalAlphabetException {
		
		List<Dirichlet> dists = new ArrayList<Dirichlet>();
		
		int minOffset = minOffset();
		int totLen = totalLength();
		int offsetTotLen = minOffset + totLen;
		
		int firstNonUniform  = totLen;
		int lastNonUniform = -1;
		boolean nonUniformFound = false;
		
		for (int i = minOffset; i < offsetTotLen; i++) {
			int pos = i - minOffset;
			Distribution[] posDists = distributionsAtPosition(i);
			if (posDists.length >= 2) {
				if (pos < firstNonUniform) {
					if (!nonUniformFound) {
						firstNonUniform = pos;
						nonUniformFound = true;
					}
				}
				lastNonUniform = pos;
				dists.add(DirichletParamEstimator
						.approximateMLE(posDists));
			} else {
				dists.add(new Dirichlet(alphabet));
			}
		}
		
		if (trimEnds) {
			for (int i = (dists.size()-1); i >= 0; i--) {
				if (i > lastNonUniform)
					dists.remove(i);
				
				if (i < firstNonUniform)
					dists.remove(i);
			}
		}
		//System.err.println("First nonuniform:" + firstNonUniform);
		//System.err.println("Last  nonuniform:" + lastNonUniform);

		MetaMotif mm = 
			new MetaMotif(dists.toArray(new Dirichlet[dists.size()]));
		mm.setName(name);
		return mm;
	}

	public Motif averageMotif() 
		throws IllegalSymbolException, IllegalAlphabetException {
		Motif motif = new Motif();
		motif.setName(name);
		
		FiniteAlphabet alphab = (FiniteAlphabet) 
			motifElems.get(0).getMotif().getWeightMatrix().getAlphabet();
		
		int minOffset = minOffset();
		int offsetTotLen = minOffset + totalLength();
		
		Distribution[] avgDists = new Distribution[totalLength()];
		
		for (int i = minOffset; i < offsetTotLen; i++) {
			Distribution avgDist = new SimpleDistribution(alphab);
			
			for (Iterator it = ((FiniteAlphabet) avgDist.getAlphabet())
					.iterator(); it.hasNext();) {
				Symbol sym = (Symbol) it.next();
				avgDist.setWeight(sym,0);
			}
			
			int motifCountAtPos = 0;
			
			for (MotifAlignmentElement me : motifElems) {
				Distribution d = me.getPosition(i);
				if (d != null) {
					motifCountAtPos++;
					for (Iterator it = ((FiniteAlphabet) 
						d.getAlphabet()).iterator(); 
						it.hasNext();) {
						Symbol sym = (Symbol) it.next();
						double w = d.getWeight(sym);
						avgDist.setWeight(
							sym, avgDist.getWeight(sym) + w);
					}
				}
			}
			
			for (Iterator 
					it = ((FiniteAlphabet) avgDist.getAlphabet())
					.iterator(); 
					it.hasNext();) {
				Symbol sym = (Symbol) it.next();
				double newW = avgDist
								.getWeight(sym)
								/(double)motifCountAtPos;
				avgDist.setWeight(sym, newW);
			}
			avgDists[i] = avgDist;
		}
		
		WeightMatrix wm = new SimpleWeightMatrix(avgDists);
		motif.setWeightMatrix(wm);
		
		return motif;
	}

	private Symbol symbolWithMaxWeight(Distribution d) {
		Symbol maxWeightSym = null;
		double maxWeight = 0;
		for (Iterator it = ((FiniteAlphabet) d.getAlphabet())
				.iterator(); it.hasNext();) {
			Symbol sym = (Symbol) it.next();
			double w;
			try {
				w = d.getWeight(sym);
			} catch (IllegalSymbolException e) {
				throw 
					new BioError(
						"Illegal symbol exception not possible here: " +
						"corrupted data structure",e);
			}
			if (w > maxWeight) {
				maxWeight = w;
				maxWeightSym = sym;
			}
		}
		return maxWeightSym;
	}
	
	private void writeOffsetAnnotations() {
		for (MotifAlignmentElement elem : motifElems) {
			Motif m = elem.getMotif();
			m.getAnnotation().setProperty("offset", elem.getOffset());
		}
	}
	
	private String consensus(Motif m) throws BioException {
		StringBuffer strBuf = new StringBuffer();
		SymbolTokenization tok = alphabet.getTokenization("token");
		
		for (int i = 0; i < m.getWeightMatrix().columns(); i++) {
			Distribution distrib = m.getWeightMatrix().getColumn(i);
			if (distrib != null) {
				strBuf.append(
						tok.tokenizeSymbol(symbolWithMaxWeight(distrib))
						);
			} else {
				strBuf.append("-");
			}
		}
		return strBuf.toString();
	}
	
	public String alignmentConsensusString() throws BioException {
		StringBuffer strBuf = new StringBuffer();
		SymbolTokenization tok = alphabet.getTokenization("token");
		
		int minOffset = minOffset();
		int maxPos = minOffset + totalLength();
		
		for (MotifAlignmentElement elem : motifElems) {
			strBuf.append(elem.getMotif().getName() + " " + elem.getOffset() + " " + "\n");
			for (int i = minOffset; i < maxPos; i++) {
				Distribution distrib = elem.getPosition(i);
				if (distrib != null) {
					strBuf.append(
							tok.tokenizeSymbol(symbolWithMaxWeight(distrib))
							);
				} else {
					strBuf.append("-");
				}
			}
			strBuf.append("\n");
		}
		
		return strBuf.toString();
	}
	
	public String consensusString() throws BioException {
		Motif avgMotif = averageMotif();
		WeightMatrix wm = avgMotif.getWeightMatrix();
		StringBuffer strBuf = new StringBuffer();
		SymbolTokenization tok = alphabet.getTokenization("token");
		
		for (int i = 0; i < wm.columns(); i++) {
			Distribution d = wm.getColumn(i);
			Symbol maxWeightSym = symbolWithMaxWeight(d);
			strBuf.append(tok.tokenizeSymbol(maxWeightSym));
		}
		
		return strBuf.toString();
	}

	public Iterator<MotifAlignmentElement> iterator() {
		return motifElems.iterator();
	}
	
	public MotifPairWithOffset highestScoringPair(MotifComparisonMatrixBundle mb)
		throws MotifAlignmentException {
		double maxScore = Double.POSITIVE_INFINITY;
		int offset = 0;
		int bestRow = -1;
		int bestCol = -1;
		
		Motif m0, m1;
		m0 = null;
		m1 = null;
		boolean flipped = false;
		
		Motif[] motifs = mb.getRowMotifs();
		for (int i = 0; i < motifs.length; i++) {
			for (int j = (i+1); j < motifs.length; j++) {
				if (alignedMotifIndices.contains(i) && 
					alignedMotifIndices.contains(j)) continue;
				//if (!mb.isAllowed(i, j)) continue;
				
				double ijScore = mb.getSenseScoreMatrix().get(i, j);
				if (ijScore < maxScore) {
					maxScore = ijScore;
					offset = mb.getSenseOffsetMatrix().get(i, j);
					bestRow = i;
					bestCol = j;
					flipped = false;
				}
				
				double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
				if (fijScore < maxScore) {
					maxScore = fijScore;
					offset = mb.getAntisenseOffsetMatrix().get(i, j);
					bestRow = i;
					bestCol = j;
					flipped = true;
				}
			}
		}
		
		if (bestRow < 0 || bestCol < 0)
			throw new MotifAlignmentException("Could not align motifs");
		
		m0 = motifs[bestRow];
		m1 = motifs[bestCol];
		
		if (m0 != null && m1 != null) {
			//System.out.println("Best row:" + bestRow + " col:" + bestCol);
			//mb.setAllowed(bestRow, bestCol, false);
			MotifPairWithOffset mpoffset = new MotifPairWithOffset(
					m0, 
					m1, 
					maxScore, 
					flipped, offset);					
			return mpoffset;
		} else {
			throw new MotifAlignmentException("Could not align motifs");
		}
	}
	
	/**
	 * Always returns the already aligned motif as the m1 of the pair 
	 * (offset is flipped accordingly)
	 * 
	 * @param mb
	 * @return
	 * @throws MotifAlignmentException
	 */
	public MotifPairWithOffset highestScoringPairWithAlreadyAligned(MotifComparisonMatrixBundle mb)
	throws MotifAlignmentException {
	//System.err.println("... comparing with already aligned");
	double maxScore = Double.POSITIVE_INFINITY;
	int offset = 0;
	int bestRow = -1;
	int bestCol = -1;
	
	Motif m0, m1;
	m0 = null;
	m1 = null;
	boolean flipped = false;
	
	Motif[] motifs = mb.getRowMotifs();
	for (int i = 0; i < motifs.length; i++) {
		for (int j = (i+1); j < motifs.length; j++) {
			/* if both are already in the alignment, let's ignore */
			if (alignedMotifIndices.contains(i) &&
					alignedMotifIndices.contains(j)) continue;
			
			/* if one of the aligned motifs is not already in the alignment, let's ignore */
			if (!alignedMotifIndices.contains(i) && 
				!alignedMotifIndices.contains(j)) continue;
			
			/*
			if (alignedMotifIndices.contains(i) &&
				alignedMotifIndices.contains(j))
				throw new IllegalStateException(
						"Both of the motifs " + allMotifs[i].getName() + 
						"(" + i + ") and " + allMotifs[j].getName() + "(" + j + ") " +
						"are already in the alignment! " +
						"This comparison should not have been made.");*/
			
			double ijScore = mb.getSenseScoreMatrix().get(i, j);
			if (ijScore < maxScore) {
				maxScore = ijScore;
				offset = mb.getSenseOffsetMatrix().get(i, j);
				bestRow = i;
				bestCol = j;
				flipped = false;
				//if (alignedMotifIndices.contains(i))
				//	iContained = true;
			}
			
			double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
			if (fijScore < maxScore) {
				maxScore = fijScore;
				offset = mb.getAntisenseOffsetMatrix().get(i, j);
				bestRow = i;
				bestCol = j;
				flipped = true;
				//if (alignedMotifIndices.contains(i))
				//	iContained = true;
			}
		}
	}
	
	if (bestRow < 0 || bestCol < 0)
		throw new MotifAlignmentException("Could not align motifs");
	
	m0 = motifs[bestRow];
	m1 = motifs[bestCol];
	
	if (m0 != null && m1 != null) {
		//System.out.println("Best row:" + bestRow + " col:" + bestCol);
		//mb.setAllowed(bestRow, bestCol, false);
		MotifPairWithOffset mpoffset;
		if (alignedMotifIndices.contains(bestRow))
			mpoffset = new MotifPairWithOffset(
				m0, 
				m1, 
				maxScore, 
				flipped, offset);
		else
			mpoffset = new MotifPairWithOffset(
				m1,
				m0,
				maxScore,
				flipped, offset);
		return mpoffset;
	} else {
		throw new MotifAlignmentException("Could not align motifs");
	}
}

	/*public MotifComparisonMatrixBundle highestScoringPair(MotifComparitorIFace mc, Motif[] motifs) 
		throws IllegalAlphabetException, IllegalSymbolException {
		return mc.getComparisonMatrixWithOffsets(motifs);
	}*/
	
	/*public MotifPairWithOffset highestScoringPair(MotifComparitorIFace mc, Motif[] motifs) 
		throws IllegalAlphabetException, IllegalSymbolException, MotifAlignmentException {
		
		MotifComparisonMatrixBundle mb = mc.
			getComparisonMatrixWithOffsets(motifs);
	
		double maxScore = Double.POSITIVE_INFINITY;
		int offset = 0;
		Motif m0, m1;
		m0 = null;
		m1 = null;
		boolean flipped = false;
		
		for (int i = 0; i < motifs.length; i++) {
			for (int j = (i+1); j < motifs.length; j++) {
				double ijScore = mb.getSenseScoreMatrix().get(i, j);
				if (ijScore < maxScore) {
					maxScore = ijScore;
					offset = mb.getSenseOffsetMatrix().get(i, j);
					m0 = motifs[i];
					m1 = motifs[j];
					flipped = false;
				}
				
				double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
				if (fijScore < maxScore) {
					maxScore = fijScore;
					offset = mb.getAntisenseOffsetMatrix().get(i, j);
					m0 = motifs[i];
					m1 = motifs[j];
					flipped = true;
				}
			}
		}
		
		if (m0 != null && m1 != null) {
			MotifPairWithOffset mpoffset = new MotifPairWithOffset(
					m0, 
					m1, 
					maxScore, 
					flipped, offset);					
			return mpoffset;
		} else {
			throw new MotifAlignmentException("Could not align motifs");
		}
	}*/
	
	/*
	public SparseMotifComparisonMatrixBundle highestScoringPairBetweenSets(
			MotifComparitorIFace mc, Motif[] motifs0, Motif[] motifs1) 
		throws Exception {
		
					mc.getComparisonMatrixWithOffsets(motifs0, motifs1));
	}*/
	
	public MotifPairWithOffset highestScoringPairBetweenSets(
			SparseMotifComparisonMatrixBundle mb) 
		throws Exception {
		Motif[] motifs0 = mb.getRowMotifs();
		Motif[] motifs1 = mb.getColumnMotifs();
			
		double maxScore = Double.POSITIVE_INFINITY;
		int offset = 0;
		Motif m0, m1;
		m0 = null;
		m1 = null;
		int bestRow, bestCol;
		bestRow = -1;
		bestCol = -1;
		
		boolean flipped = false;
		
		//System.err.println("---");
		for (int i = 0; i < motifs0.length; i++) {
			for (int j = 0; j < motifs1.length; j++) {
				if (!mb.isAllowed(i, j)) continue;
				
				double ijScore = mb.getSenseScoreMatrix().get(i, j);
				if (ijScore < maxScore) {
					maxScore = ijScore;
					offset = mb.getSenseOffsetMatrix().get(i, j);
					//System.err.println("score::"+ijScore + " offset:" + offset);
					bestRow = i;
					bestCol = j;
					flipped = false;
				}
				
				double fijScore = mb.getAntisenseScoreMatrix().get(i, j);
				if (fijScore < maxScore) {
					maxScore = fijScore;
					offset = mb.getAntisenseOffsetMatrix().get(i, j);
					//System.err.println("f score::"+ijScore + " offset:" + offset);
					m0 = motifs0[i];
					m1 = motifs1[j];
					bestRow = i;
					bestCol = j;
					
					flipped = true;
				}
			}
		}
		
		m0 = motifs0[bestRow];
		m1 = motifs1[bestCol];
		
		
		if (m0 != null && m1 != null) {
			mb.setAllowed(bestRow, bestCol, false);
			MotifPairWithOffset mpoffset = 
				new MotifPairWithOffset(m0, m1, maxScore, flipped, offset);					
			return mpoffset;
		} else {
			throw new MotifAlignmentException("Could not align motifs");
		}

	}

	
	private int firstPosWithMinNumCols(int minCols) throws NoPositionFoundException {
		int minOffset = minOffset();
		int maxPos = minOffset + totalLength();
		
		for (int i = minOffset; i < maxPos; i++)
			if (distributionsAtPosition(i).length >= minCols) return i;
		
		throw new NoPositionFoundException(
				"No positions with the minimum number of " + 
				minCols + " columns defined");
	}
	
	private int lastPosWithMinNumCols(int minCols) throws NoPositionFoundException {
		int minOffset = minOffset();
		int maxPos = minOffset + totalLength();
		
		for (int i = maxPos; i > minOffset; i--)
			if (distributionsAtPosition(i).length >= minCols) return i;

		throw new NoPositionFoundException(
				"No positions with the minimum number of " + 
				minCols + " columns defined");
	}
	
	private int firstPosWithMinNumColsForElem(MotifAlignmentElement elem) {
		int minOffset = minOffset();
		int maxPos = minOffset + totalLength();
		
		for (int i = minOffset; i < maxPos; i++)
			if (elem.getPosition(i) != null)
				return i;
		
		throw new MotifAlignmentError(
			"The element " + elem + " " +
			"has no columns defined in any of the alignment positions. " +
			"The alignment data structure is corrupted.");
		
	}
	
	public MotifAlignment fillMissingColumns(Distribution fillerColumn) 
		throws IllegalAlphabetException {
		MotifAlignment alignment = new MotifAlignment(alphabet);
		
		int minOffset = minOffset();
		int maxPos = minOffset + totalLength();
		int totLen = totalLength();
		
		for (MotifAlignmentElement elem : motifElems) {
			Distribution[] dists = new Distribution[totLen];
			for (int i = minOffset; i < maxPos; i++) {
				Distribution d = elem.getPosition(i);
				if (d == null)
					dists[i] = new SimpleDistribution(fillerColumn);
				else
					dists[i] = new SimpleDistribution(d);
			}
			WeightMatrix wm = new SimpleWeightMatrix(dists);
			Motif m = new Motif(elem.getMotif());
			m.setWeightMatrix(wm);
			alignment.addMotif(m, 0, elem.isFlipped());
		}
		
		return alignment;
	}
	
	public MotifAlignment trimToColumnsPerPosition(int minColPerPos) throws NoPositionFoundException {
		MotifAlignment alignment = new MotifAlignment(alphabet);
		
		if (minColPerPos < 1) {
			throw new IllegalArgumentException("Allowed argument range: minColPerPos >= 1");
		}
		if (minColPerPos == 1) {
			return new MotifAlignment(this);
		}
		
		int first = firstPosWithMinNumCols(minColPerPos);
		int last = lastPosWithMinNumCols(minColPerPos);
		int minOffset = minOffset();
		
		for (MotifAlignmentElement elem : motifElems) {
			List<Distribution> dists = new ArrayList<Distribution>();
			boolean firstOneFound = false;
			int addedOffset = 0;
			for (int i = first; i <= last; i++) {
				Distribution d = elem.getPosition(i);
				if (d != null) {
					if (!firstOneFound) firstOneFound = true;
					dists.add(d);
				} else if (!firstOneFound) {
					addedOffset++;
				}
			}
			
			if (dists.size() > 0) {
				WeightMatrix wm;
				try {
					wm = new SimpleWeightMatrix(
							dists.toArray(
								new Distribution[dists.size()]));
				} catch (IllegalAlphabetException e) {
					throw new BioError(
							"Illegal alphabet exception should not occur here. " +
							"Corrupted data structure!",e);
				}
				Motif m = new Motif(elem.getMotif());
				m.setWeightMatrix(wm);
				alignment.addMotif(m, 
									addedOffset, 
									elem.isFlipped());
			}
		}
		
		return alignment;
	}
}
