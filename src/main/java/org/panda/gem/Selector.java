package org.panda.gem;

import org.panda.utility.statistics.FDR;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Selects significant and categorized triplets.
 *
 * @author Ozgun Babur
 */
public class Selector
{
	/**
	 * Selects significant and categorized triplets.
	 * @param trips triplets
	 * @param fdrThr FDR threshold for selecting gamma and then betaM
	 * @param categThr P-value threshold to use during category assignment
	 * @return list of significant and categorized triplets
	 */
	public static List<Triplet> selectSignificantAndCategorized(Collection<Triplet> trips,
		double fdrThr, double categThr)
	{
		// select with gamma pval

		Map<Triplet, Double> map = new HashMap<>(trips.size());

		trips.stream()
			.peek(Triplet::initGamma)
			.forEach(t -> map.put(t, t.gamma.p));

		List<Triplet> tripList = FDR.selectBH(map, fdrThr);

		// select with betaM

		map.clear();

		tripList.stream()
			.peek(Triplet::initBetaM)
			.forEach(t -> map.put(t, t.betaM.p));

		tripList = FDR.selectBH(map, fdrThr);

		// select with alphaM/betaM, and being in a category

		tripList = tripList.stream()
			.peek(Triplet::initOtherCoefficients)
			.filter(t -> t.alphaM.v / t.betaM.v < 1)
			.peek(t -> t.initCategory(categThr))
			.filter(t -> t.cat != null)
			.collect(Collectors.toList());

		return tripList;
	}
}
