package org.panda.gem.resource;

import org.panda.gem.Gene;
import org.panda.gem.ModPrint;
import org.panda.gem.Selector;
import org.panda.gem.Triplet;
import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.Summary;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * This class is designed for using GEO files in GEM analysis. It specifically handles the expO dataset in GEO. The user
 * should manually download the related platform file (GPL...) and the series file (GSE...) from GEO using its web
 * interface. Then use pass those files to the constructor of this class. The main method of this class is also an
 * example for its usage.
 *
 * @author Ozgun Babur
 */
public class GEOFromFileLoader implements GeneProvider
{
	/**
	 * Name of the platform file.
	 */
	String platformFile;

	/**
	 * Name of the series file.
	 */
	String valuesFile;

	/**
	 * Cache for not creating redundant genes.
	 */
	private Map<String, Gene> cache;

	/**
	 * Map from row ID to the value array.
	 */
	private Map<String, double[]> id2vals;

	/**
	 * Map from gene symbol to related row IDs.
	 */
	private Map<String, Set<String>> sym2IDs;

	/**
	 * Constructor with two necessary GEO files.
	 */
	public GEOFromFileLoader(String platformFile, String valuesFile) throws IOException
	{
		this.platformFile = platformFile;
		this.valuesFile = valuesFile;
		readPlatform();
		readValues();
		cache = new HashMap<>();
	}

	/**
	 * Reads the platform file and prepares the mapping of symbols to related row IDs.
	 */
	private void readPlatform() throws IOException
	{
		sym2IDs = new HashMap<>();

		String[] header = Files.lines(Paths.get(platformFile)).filter(l -> l.startsWith("ID\t")).map(l -> l.split("\t"))
			.findFirst().get();

		int symbolIndex = ArrayUtil.indexOf(header, "Gene Symbol");

		Files.lines(Paths.get(platformFile)).filter(l -> !l.isEmpty()).filter(l -> !l.startsWith("!"))
			.filter(l -> !l.startsWith("^")).filter(l -> !l.startsWith("ID\t")).map(l -> l.split("\t"))
			.filter(t -> t.length > symbolIndex).filter(t -> !t[symbolIndex].isEmpty()).forEach(t ->
		{
			String id = t[0];
			String[] syms = t[symbolIndex].split(" /// ");
			for (String sym : syms)
			{
				if (!sym2IDs.containsKey(sym)) sym2IDs.put(sym, new HashSet<>());
				sym2IDs.get(sym).add(id);
			}
		});
	}

	/**
	 * Reads the values file and prepared the id to value array mapping.
	 */
	private void readValues() throws IOException
	{
		id2vals = new HashMap<>();

		Files.lines(Paths.get(valuesFile)).filter(l -> !l.startsWith("!")).filter(l -> !l.isEmpty())
			.filter(l -> !l.startsWith("\"ID_REF\"\t")).map(l -> l.split("\t")).forEach(t ->
		{
			double[] v = new double[t.length - 1];
			String id = t[0].replaceAll("\"", "");
			for (int i = 1; i < t.length; i++)
			{
				v[i - 1] = Double.valueOf(t[i]);
			}
			id2vals.put(id, v);
		});
	}

	@Override
	public Gene get(String symbol)
	{
		if (cache.containsKey(symbol)) return cache.get(symbol);

		if (sym2IDs.containsKey(symbol))
		{
			List<double[]> valList = new ArrayList<>();
			for (String id : sym2IDs.get(symbol))
			{
				valList.add(id2vals.get(id));
			}

			Gene gene = new Gene(symbol, selectHighestVariation(valList));
			cache.put(symbol, gene);
			return gene;
		}

		cache.put(symbol, null);
		return null;
	}

	/**
	 * For each gene symbol, the row with highest variance is selected. This method selects that row.
	 */
	private double[] selectHighestVariation(List<double[]> list)
	{
		if (list.size() == 1) return list.get(0);

		double maxVar = 0;
		double[] max = null;

		for (double[] vals : list)
		{
			double var = Summary.variance(vals);

			if (var > maxVar)
			{
				maxVar = var;
				max = vals;
			}
		}

		return max;
	}
}
