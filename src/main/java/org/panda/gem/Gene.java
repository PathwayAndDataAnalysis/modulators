package org.panda.gem;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represents a gene with an HGNC symbol, and an expression array. It also has a discretized expression where each 3
 * discrete value represents a tertile.
 *
 * @author Ozgun Babur
 */
public class Gene
{
	/**
	 * HGNC symbol.
	 */
	public String symbol;

	/**
	 * Expression array.
	 */
	public double[] vals;

	/**
	 * Discretized expression array.
	 * 0: low 1: high -1: middle
	 */
	int[] status;

	public Gene(String symbol, double[] vals)
	{
		this.symbol = symbol;

		if (vals != null)
		{
			this.vals = vals;
			assignStatus();
		}
	}

	/**
	 * Discretizes values according to corresponding tertiles.
	 */
	private void assignStatus()
	{
		status = new int[vals.length];

		// indices sorted to values
		List<Integer> indices = IntStream.range(0, vals.length).boxed()
			.sorted((i1, i2) -> (int) Math.signum(vals[i1] - vals[i2])).collect(Collectors.toList());

		int m1 = vals.length / 3;
		int m2 = (vals.length * 2) / 3;

		for (int i = 0; i < m1; i++) status[indices.get(i)] = 0;
		for (int i = m1; i < m2; i++) status[indices.get(i)] = -1;
		for (int i = m2; i < vals.length; i++) status[indices.get(i)] = 1;
	}
}
