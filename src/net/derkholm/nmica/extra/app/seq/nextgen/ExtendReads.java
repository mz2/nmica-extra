package net.derkholm.nmica.extra.app.seq.nextgen;

import java.util.ArrayList;
import java.util.List;

import net.derkholm.nmica.build.NMExtraApp;
import net.derkholm.nmica.build.VirtualMachine;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import org.biojava.bio.BioException;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@NMExtraApp(launchName = "ngextend", vm = VirtualMachine.SERVER)
@App(overview = "Extend reads by specified number of nucleotides. " +
		"Note that the output simply moves the read start along the reference " +
		"and adds elements to the cigar string. " +
		"That is, the sequence read itself is not changed and each of the SAM records is in fact " +
		"invalid.", generateStub = true)
public class ExtendReads extends FilteringSAMProcessor {
	private int extraCigarLength;
	
	@Option(help="Expand reads by specified number of nucleotides (bound by reference sequence ends)")
	public void setBy(int i) {
		this.extraCigarLength = i;
	}
	
	public void main(String[] args) throws BioException {
		setIterationType(IterationType.ONE_BY_ONE);
		initializeSAMReader();
		initializeSAMWriter(this.sorted);
		
		process();
		outWriter.close();
	}
	
	@Override
	public void process(SAMRecord rec, int readIndex) {
		int len = refSeqLengths.get(rec.getReferenceName());
		Cigar oldCigar = rec.getCigar();
		List<CigarElement> cigarElems = new ArrayList<CigarElement>(oldCigar.getCigarElements());
		
		/* positive strand */
		if (!rec.getReadNegativeStrandFlag()) {
			int extraLen = extraCigarLength;
			int distFromRefEnd = len - (rec.getAlignmentStart() + rec.getCigarLength());
			
			if (distFromRefEnd < extraCigarLength) {
				extraLen = distFromRefEnd - 1;
			}
			
			cigarElems.add(new CigarElement(extraLen,CigarOperator.M));
			rec.setCigar(new Cigar(cigarElems));
			
			//rec.setAlignmentEnd(Math.min(len, rec.getAlignmentEnd() + extraCigarLength));
		} 
		/* negative strand */
		else {
			int extraLen = extraCigarLength;
			int distFromRefStart = rec.getAlignmentStart();
			
			if (distFromRefStart - extraCigarLength < 0) {
				extraLen = distFromRefStart;
			}
			
			cigarElems.add(0,new CigarElement(extraLen,CigarOperator.M));
			rec.setCigar(new Cigar(cigarElems));
			
			rec.setAlignmentStart(Math.max(0, distFromRefStart - extraCigarLength));
		}
		
		outWriter.addAlignment(rec);
	}
}