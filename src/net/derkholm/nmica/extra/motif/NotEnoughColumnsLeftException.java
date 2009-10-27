package net.derkholm.nmica.extra.motif;

import org.biojava.bio.BioException;

public class NotEnoughColumnsLeftException extends BioException {
	private static final long serialVersionUID = -4359774802811253827L;

	public NotEnoughColumnsLeftException(String str) {
		super(str);
	}
	
	public NotEnoughColumnsLeftException() {
		super();
	}
}