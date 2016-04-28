package org.panda.gem;

import org.panda.utility.Tuple;

/**
 * Defines and detects 6 modulation categories.
 * @author Ozgun Babur
 */
public enum ModulationCategory
{
	INVERTS_ACTIVATION(-1, 1, -1, -1, -1),
	INVERTS_INHIBITION(1, -1, 1, 1, 1),
	ENHANCES_ACTIVATION(1, 0, 1, 1, 1),
	ENHANCES_INHIBITION(-1, 0, -1, -1, -1),
	ATTENUATES_ACTIVATION(-1, 1, 0, 0, 0),
	ATTENUATES_INHIBITION(1, -1, 0, 0, 0);

	/**
	 * Requirement vector indicating signs of the triplet coefficients.
	 * @see Triplet#getCoefficients()
	 */
	int[] r;

	ModulationCategory(int... r)
	{
		this.r = r;
	}

	/**
	 * Finds the category that matches the coefficient array of a triplet. Note that each invert category also matches
	 * to two other categories, so inversions are checked first, and others are matched only if inversions did not
	 * match.
	 * @param tuples the coefficient array
	 * @param thr p-value threshold for a significant coefficient
	 * @return the matching category
	 */
	private boolean matches(Tuple[] tuples, double thr)
	{
		assert tuples.length == r.length;

		for (int i = 0; i < r.length; i++)
		{
			if (r[i] != 0)
			{
				if (tuples[i].p > thr) return false;
				if (Math.signum(tuples[i].v) * r[i] < 0) return false;
			}
		}
		return true;
	}

	/**
	 * Checks if this category matches the triplet coefficients.
	 */
	public static ModulationCategory match(Triplet t, double thr)
	{
		for (ModulationCategory cat : values())
		{
			if (cat.matches(t.getCoefficients(), thr)) return cat;
		}
		return null;
	}
}
