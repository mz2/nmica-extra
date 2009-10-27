package biobits.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

import biobits.utils.Collects;
import biobits.utils.Function;

public class IOTools {
	private IOTools() {
	}
	
	public static BufferedReader fileBufferedReader(File f) 
		throws Exception
	{
		return new BufferedReader(fileReader(f));
	}
	
	public static Reader fileReader(File f)
		throws Exception
	{
		if (f.getName().endsWith(".gz")) {
			return new InputStreamReader(new GZIPInputStream(new FileInputStream(f)));
		} else {
			return new FileReader(f);
		}
	}

	private static Reader inputReader(String arg)
		throws Exception
	{
		if ("-".equals(arg)) {
			return new InputStreamReader(System.in);
		} else {
			return fileReader(new File(arg));
		}
	}
	
	public static Reader inputReader(String[] args) 
		throws Exception
	{
		if (args.length == 0) {
			return inputReader("-");
		} else if (args.length == 1) {
			return inputReader(args[0]);
		} else {
			return new SequenceReader(
					Collects.map(
							new Function<String,Reader>() {
									public Reader apply(String name) {
										try {
											return inputReader(name);
										} catch (Exception ex) {
											throw new RuntimeException(ex);
										}
									}
								
							},
							Arrays.asList(args)
					)
			);
		}
	}
	
	public static BufferedReader inputBufferedReader(String[] args) 
		throws Exception
	{
		return new BufferedReader(inputReader(args));
	}
}
