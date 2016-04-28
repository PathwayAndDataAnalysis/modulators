package org.panda.gem;

import org.panda.gem.resource.GeneProvider;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.ErrorFunction;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Modulator - Factor - Target triplet of genes. This class holds the case counts and other statistics about the
 * triplet.
 *
 * @author Ozgun Babur
 */
public class Triplet
{
	public static final double SQRT2 = Math.sqrt(2);

	/**
	 * Modulator
	 */
	public Gene M;

	/**
	 * Transcription Factor
	 */
	public Gene F;

	/**
	 * Target
	 */
	public Gene T;

	/**
	 * Frequencies of 8 corner cases.
	 */
	int[][][] f;

	/**
	 * Totals for 4 M-F status.
	 */
	int[][] n;

	/**
	 * Proportions for 4 M-F status.
	 */
	double[][] p;

	// Inferred coefficients of the triplet equation, along with their significances

	Tuple alphaM;
	Tuple betaM;
	Tuple alphaF;
	Tuple betaF;
	Tuple gamma;

	/**
	 * AlphaF + BetaM
	 */
	Tuple aFbM;

	/**
	 * The detected modulation category.
	 */
	public ModulationCategory cat;

	public Triplet(Gene m, Gene f, Gene t)
	{
		M = m;
		F = f;
		T = t;

		if (hasAllGenes())
		{
			initFreqs();
			initTotalsAndProportions();
		}
	}

	public boolean hasAllGenes()
	{
		return M != null && F != null && T != null;
	}

	public Triplet(String lineOfFile, GeneProvider loader)
	{
		String[] token = lineOfFile.split("\t");
		M = loader.get(token[0]);
		F = loader.get(token[1]);
		T = loader.get(token[2]);
		cat = ModulationCategory.valueOf(token[3]);

		f = new int[2][2][2];

		int x = 4;

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					f[i][j][k] = Integer.parseInt(token[x++]);

		initTotalsAndProportions();
	}

	/**
	 * Gets the coefficients in an array for the modulation categories to evaluate their match.
	 */
	public Tuple[] getCoefficients()
	{
		return new Tuple[]{gamma, alphaF, betaF, betaM, aFbM};
	}

	private void initFreqs()
	{
		f = new int[2][2][2];
		for (int i = 0; i < M.status.length; i++)
		{
			if (M.status[i] < 0 || F.status[i] < 0 || T.status[i] < 0)
				continue;
			f[M.status[i]][F.status[i]][T.status[i]]++;
		}
	}

	public void initTotalsAndProportions()
	{
		n = new int[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				n[i][j] = f[i][j][0] + f[i][j][1];

		p = new double[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				p[i][j] = f[i][j][1] / (double) n[i][j];
	}

	public void initGamma()
	{
		double pg = (f[0][0][1] + f[0][1][1] + f[1][0][1] + f[1][1][1]) /
			(double) (n[0][0] + n[0][1] + n[1][0] + n[1][1]);

		double stdev = Math.sqrt(pg * (1 - pg) * ((1D / n[0][0]) + (1D / n[0][1]) + (1D / n[1][0]) + (1D / n[1][1])));

		gamma = generateTuple(p[1][1] - p[0][1] - p[1][0] + p[0][0], stdev);
	}

	public void initBetaM()
	{
		betaM = pairwise(1, 1, 0, 1);
	}

	public void initOtherCoefficients()
	{
		alphaM = pairwise(1, 0, 0, 0);
		alphaF = pairwise(0, 1, 0, 0);
		betaF = pairwise(1, 1, 1, 0);
		aFbM = pairwise(1, 1, 0, 0);
	}

	public void initCategory(double thr)
	{
		cat = ModulationCategory.match(this, thr);
	}

	private Tuple pairwise(int i, int j, int k, int l)
	{
		double v = p[i][j] - p[k][l];
		double pijkl = (f[i][j][1] + f[k][l][1]) / (double) (n[i][j] + n[k][l]);
		double stdev = Math.sqrt(pijkl * (1 - pijkl) * ((1D / n[i][j]) + (1D / n[k][l])));
		return generateTuple(v, stdev);
	}

	private Tuple generateTuple(double v, double stdev)
	{
		return new Tuple(v, 1 - ErrorFunction.getSignif(Math.abs(v) / (SQRT2 * stdev)));
	}

	// Section: Static File Operations

	/**
	 * Loads list of triplets from a file.
	 */
	public static List<Triplet> load(String file) throws IOException
	{
		DummyGeneProvider loader = new DummyGeneProvider();
		return Files.lines(Paths.get(file)).skip(1).map(line -> new Triplet(line, loader)).collect(Collectors.toList());
	}

	/**
	 * Column headers of the tab-delimited triplet file.
	 */
	public static final String FILE_HEADER =
		"Modulator\tFactor\tTarget\tModulation Category\t000\t001\t010\t011\t100\t101\t110\t111";

	/**
	 * Writes list of triplets to a file.
	 */
	public static void write(List<Triplet> trips, String file) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		writer.write(FILE_HEADER);
		trips.forEach(t -> FileUtil.write("\n" + t, writer));
		writer.close();
	}


	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder(ArrayUtil.getString("\t", M.symbol, F.symbol, T.symbol, cat));
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					sb.append("\t").append(f[i][j][k]);
		return sb.toString();
	}

	private String key()
	{
		return M.symbol + F.symbol + T.symbol;
	}

	public static class DummyGeneProvider implements GeneProvider
	{
		private Map<String, Gene> cache = new HashMap<>();

		@Override
		public Gene get(String symbol)
		{
			if (!cache.containsKey(symbol)) cache.put(symbol, new Gene(symbol, null));
			return cache.get(symbol);
		}
	}

	@Override
	public int hashCode()
	{
		return M.symbol.hashCode() + F.symbol.hashCode() + T.symbol.hashCode() + (cat == null ? 0 : cat.hashCode());
	}

	@Override
	public boolean equals(Object o)
	{
		return o instanceof Triplet &&
			((Triplet) o).M.symbol.equals(M.symbol) &&
			((Triplet) o).F.symbol.equals(F.symbol) &&
			((Triplet) o).T.symbol.equals(T.symbol) &&
			((Triplet) o).cat == cat;
	}
}
