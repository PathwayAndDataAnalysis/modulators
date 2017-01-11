package org.panda.gem.resource;

import org.panda.gem.Gene;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * This class is for loading a simple expression file, where the first row is header, first column is the gene symbol,
 * and other cells contain the expression values.
 *
 * @author Ozgun Babur
 */
public class SimpleFileExpressionLoader implements GeneProvider
{
	/**
	 * Name of the tab-delimited expression file.
	 */
	String filename;

	/**
	 * Cache for not creating redundant genes.
	 */
	private Map<String, Gene> cache;

	/**
	 * Constructor with the filename.
	 */
	public SimpleFileExpressionLoader(String filename) throws IOException
	{
		this.filename = filename;
		cache = new HashMap<>();
		readFile();
	}

	/**
	 * Reads the values file and load genes in a cache.
	 */
	private void readFile() throws IOException
	{
		Files.lines(Paths.get(filename)).filter(l -> !l.startsWith("!")).filter(l -> !l.isEmpty())
			.filter(l -> !l.startsWith("#")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			double[] v = new double[t.length - 1];
			String symbol = t[0].replaceAll("\"", "");
			if (symbol.contains("|")) symbol = symbol.substring(0, symbol.indexOf("|"));
			for (int i = 1; i < t.length; i++)
			{
				v[i - 1] = Double.valueOf(t[i]);
			}
			cache.put(symbol, new Gene(symbol, v));
		});
	}

	@Override
	public Gene get(String symbol)
	{
		return cache.get(symbol);
	}
}
