package net.derkholm.nmica.motif.align;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import net.derkholm.nmica.motif.Motif;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

public class MotifAlignment {
	protected List<MotifAlignmentElement> motifs;
	protected SortedSet<MotifAlignmentElement> offsetSortedMotifs;
	protected HashMap<Motif, MotifAlignmentElement> motifsToAlignmentElems;
	
	protected String name;
	protected final FiniteAlphabet alphabet;
	
	
	public MotifAlignment(FiniteAlphabet alphabet) {
		this.alphabet = alphabet;
		
		motifs = new ArrayList<MotifAlignmentElement>();
		offsetSortedMotifs = 
			new TreeSet<MotifAlignmentElement>(
					new MotifAlignmentElementOffsetComparator());
		
		motifsToAlignmentElems = 
			new HashMap<Motif, MotifAlignmentElement>();
	}
		
	public void addMotif(Motif m, int offset) {
		System.out.println("+"+m.getName()+" offset:" + offset);
		if (!contains(m)) {
			MotifAlignmentElement me = new MotifAlignmentElement(m,offset);
			motifs.add(me);
			motifsToAlignmentElems.put(m, me);
			offsetSortedMotifs.add(me);
		} else {
			throw new IllegalArgumentException("Motif "+m.getName()+" is already present in the alignment (offset " + offset + ")");
		}
	}
	
	public boolean removeMotif(Motif m) {
		if (motifsToAlignmentElems.containsKey(m)) {
			MotifAlignmentElement me = motifsToAlignmentElems.get(m);
			motifs.remove(me);
			offsetSortedMotifs.remove(me);
			motifsToAlignmentElems.remove(m);
			return true;
		} else {
			return false;
		}
	}
	
	public boolean contains(Motif m) {
		return motifs.contains(motifsToAlignmentElems.get(m));
	}
	
	public int offset(Motif m) {
		return motifsToAlignmentElems.get(m).getOffset();
	}

	public Iterator<MotifAlignmentElement> iterator() {
		return motifs.iterator();
	}

	public int numberOfMotifs() {
		return motifs.size();
	}
	
	public int totalLength() {
		MotifAlignmentElement mel = offsetSortedMotifs.last();
		MotifAlignmentElement mef = offsetSortedMotifs.first();
		
		return mel.getOffset() + mel.getMotif().getWeightMatrix().columns() - mef.getOffset();
	}

	public MotifAlignmentElement[] toArray(MotifAlignmentElement[] motifs) {
		return this.motifs.toArray(motifs);
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Motif averageMotif() 
		throws IllegalSymbolException, IllegalAlphabetException {
		Motif motif = new Motif();
		motif.setName(name);
		
		FiniteAlphabet alphab = (FiniteAlphabet) 
			motifs.get(0).getMotif().getWeightMatrix().getAlphabet();
		
		int minOffset = minOffset();
		int maxOffset = maxOffset();
		Distribution[] avgDists = new Distribution[totalLength()];
		
		for (int i = minOffset; i < maxOffset; i++) {
			Distribution avgDist = new SimpleDistribution(alphab);
			double sum = 0;
			
			for (int m = 0; i < motifs.size(); m++) {
				MotifAlignmentElement me = motifs.get(i);
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

	private int maxOffset() {
		return offsetSortedMotifs.last().getOffset();
	}

	private int minOffset() {
		return offsetSortedMotifs.first().getOffset();
	}

	public class MotifAlignmentElement {
		private Motif motif;
		private int offset;
		
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

		public MotifAlignmentElement(Motif m, int offset) {
			this.motif = m;
			this.offset = offset;
		}
		
		public Distribution getPosition(int position) {
			int offsetPos = offset + position;
			
			if (motif.getWeightMatrix().columns() < offsetPos && offsetPos >= 0) {
				return motif.getWeightMatrix().getColumn(offsetPos);
			} else
				return null;
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
		Motif[] ms = new Motif[motifs.size()];
		for (int i = 0; i < motifs.size(); i++)
			ms[i] = motifs.get(i).getMotif();
		
		return ms;
	}

}
