package biobits.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

public class Text {
	private Text() {}
	
	public static List<String> accumulate(StringTokenizer toke) {
		List<String> l = new ArrayList<String>();
		while (toke.hasMoreTokens()) {
			l.add(toke.nextToken());
		}
		return l;
	}
	
	public static List<String> split(String s) {
		return accumulate(new StringTokenizer(s));
	}
	
	public static List<String> split(String s, String tokChars) {
		return accumulate(new StringTokenizer(s, tokChars));
	}

	public static String remove(String s, String killChars) {
		StringBuilder sb = new StringBuilder();
		for (int c = 0; c < s.length(); ++c) {
			char cc = s.charAt(c);
			if (killChars.indexOf(cc) < 0) {
				sb.append(cc);
			}
		}
		return sb.toString();
	}
}
