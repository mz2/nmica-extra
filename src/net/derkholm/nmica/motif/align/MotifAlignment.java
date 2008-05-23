package net.derkholm.nmica.motif.align;

import java.beans.PropertyChangeSupport;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import net.derkholm.nmica.model.metamotif.DirichletParamEstimator;
import net.derkholm.nmica.model.metamotif.MetaMotif;
import net.derkholm.nmica.motif.Motif;

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
	
	private final PropertyChangeSupport pcs = new PropertyChangeSupport(this);
	
	protected String name;
	protected final FiniteAlphabet alphabet;
	
	
	public MotifAlignment(FiniteAlphabet alphabet) {
		this.alphabet = alphabet;
		
		motifElems = new ArrayList<MotifAlignmentElement>();
		offsetSortedMotifs = 
			new TreeSet<MotifAlignmentElement>(
					new MotifAlignmentElementOffsetComparator());
		
		motifsToAlignmentElems = 
			new HashMap<Motif, MotifAlignmentElement>();
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
		//System.err.println("+"+m.getName()+" offset:" + offset);
		if (!contains(m)) {
			MotifAlignmentElement me = new MotifAlignmentElement(m,offset,flipped);
			motifElems.add(me);
			motifsToAlignmentElems.put(m, me);
			offsetSortedMotifs.add(me);
		} else {
			throw new IllegalArgumentException("Motif "+m.getName()+" is already present in the alignment (offset " + offset + ")");
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
		return motifsToAlignmentElems.get(m).getOffset();
	}

	public int numberOfMotifs() {
		return motifElems.size();
	}
	
	public int totalLength() {
		MotifAlignmentElement mel = offsetSortedMotifs.last();
		MotifAlignmentElement mef = offsetSortedMotifs.first();
		
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
		
		/*
		return mel.getOffset() + 
				mel.getMotif().getWeightMatrix().columns() - 
				mef.getOffset();*/
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

	public Motif[] motifs() {
		writeOffsetAnnotations();
		Motif[] ms = new Motif[motifElems.size()];
		for (int i = 0; i < motifElems.size(); i++)
			ms[i] = motifElems.get(i).getMotif();
		
		return ms;
	}

	
	public Distribution[] distributionsAtPosition(int position) {
		Distribution[] dists = new Distribution[this.motifElems.size()];
		
		int i = 0;
		for (MotifAlignmentElement elem : motifElems) {
			Distribution d = elem.getPosition(position);
			assert d != null;
			dists[i++] = d;
		}
		
		return dists;
	}
	
	public MetaMotif metamotif() throws IllegalSymbolException, IllegalAlphabetException {
		Distribution[] dists = new Distribution[totalLength()];
		for (int i = minOffset(); i <= maxOffset(); i++) {
			dists[i] = DirichletParamEstimator.approximateMLE(distributionsAtPosition(i));
		}
		
		MetaMotif mm = new MetaMotif(dists);
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
		int maxOffset = maxOffset();
		Distribution[] avgDists = new Distribution[totalLength()];
		
		for (int i = minOffset; i < maxOffset; i++) {
			Distribution avgDist = new SimpleDistribution(alphab);
			double sum = 0;
			
			for (int m = 0; i < motifElems.size(); m++) {
				MotifAlignmentElement me = motifElems.get(i);
				Distribution d = me.getPosition(m);
					
				for (Iterator it = ((FiniteAlphabet) d.getAlphabet())
						.iterator(); it.hasNext();) {
					Symbol sym = (Symbol) it.next();
					double w = d.getWeight(sym);
					avgDist.setWeight(sym, avgDist.getWeight(sym) + w);
					sum = sum + w;
				}
			}
			
			for (Iterator it = ((FiniteAlphabet) avgDist.getAlphabet())
					.iterator(); it.hasNext();) {
				Symbol sym = (Symbol) it.next();
				avgDist.setWeight(sym, avgDist.getWeight(sym)/sum);
			}
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
	
	public String alignmentConsensusString() throws BioException {
		StringBuffer strBuf = new StringBuffer();
		SymbolTokenization tok = alphabet.getTokenization("token");
		
		System.err.println("Min offset:" + minOffset());
		System.err.println("Max offset:" + maxOffset());
		System.err.println("Total length:" + totalLength());
		
		int minOffset = minOffset();
		int maxPos = minOffset + totalLength();
		
		for (MotifAlignmentElement elem : motifElems) {
			strBuf.append(elem.getMotif().getName() + "\n");
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

}
