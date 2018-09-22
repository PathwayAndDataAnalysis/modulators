package org.panda.gem.resource;

import org.panda.gem.Gene;
import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.UniquePrinter;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Loads expression from a downloaded TCGA dataset.
 * @author Ozgun Babur
 */
public class TCGAExpressionLoader implements GeneProvider
{
	/**
	 * Cache for not creating redundant genes.
	 */
	private Map<String, Gene> cache;

	/**
	 * TCGA expression reader.
	 */
	private ExpressionReader expR;

	/**
	 * Available samples in the TCGA dataset.
	 */
	private String[] samples;

	private double stdevThr;

	private UniquePrinter up = new UniquePrinter();

	public TCGAExpressionLoader(String dirForExpressions, Set<String> sampleSubset) throws FileNotFoundException
	{
		expR = new ExpressionReader(dirForExpressions + "/expression.txt", null, 15);
		samples = sampleSubset == null ? expR.getSamples().stream().sorted().toArray(String[]::new) :
			expR.getSamples().stream().filter(sampleSubset::contains).sorted().toArray(String[]::new);
		cache = new HashMap<>();
		stdevThr = 0;
	}

	@Override
	public Gene get(String symbol)
	{
		if (!cache.containsKey(symbol))
		{
			double[] vals = expR.getGeneAlterationArray(symbol, samples);
			if (vals != null)
			{
				if (Summary.stdev(vals) >= stdevThr)
				{
					Gene g = new Gene(symbol, vals);
					cache.put(symbol, g);
				}
				else
				{
					cache.put(symbol, null);
					up.print("Gene does not pass stdev threshold = ", symbol);
				}
			}
			else
			{
				cache.put(symbol, null);
				up.print("Gene is not found = ", symbol);
			}
		}
		return cache.get(symbol);
	}

	public void setStdevThr(double stdevThr)
	{
		this.stdevThr = stdevThr;
	}

	public void writeExpressionHistograms()
	{
		try
		{
			OutputStream os = new FileOutputStream("/home/babur/Documents/GEM/TCGAPC/temp.txt");

			for (String gene : cache.keySet())
			{
				Histogram h = new Histogram(1, gene);
				h.setBorderAtZero(true);
				Gene g = cache.get(gene);
				if (g != null)
				{
					h.countAll(g.vals);
					h.write(os);
				}
			}
			os.close();

		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}
}
