package net.derkholm.nmica.motif.align;

public class MotifAlignmentException extends Exception {
	public MotifAlignmentException(String string) {
		super(string);
	}
	
	public MotifAlignmentException(String s, Exception e) {
		super(s,e);
	}
	
	public MotifAlignmentException(Exception e) {
		super(e);
	}

	private static final long serialVersionUID = 1638646538172596500L;

}
