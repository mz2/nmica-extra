package net.derkholm.nmica.motif.align;

public class MotifAlignmentError extends Error {
	private static final long serialVersionUID = -4042625236433483460L;

	public MotifAlignmentError() {super();}
	public MotifAlignmentError(String str) {super(str);}
	public MotifAlignmentError(String str, Exception e) {super(str, e);}
}
