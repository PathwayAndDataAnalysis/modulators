package org.panda.gem.resource;

import org.panda.gem.Triplet;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Generates triplets based on the user-provided factor, interacting proteins, and targets.
 *
 * @author Ozgun Babur
 */
public class CustomTripletMaker
{
	/**
	 * Generates custom triplets for the given factor.
	 */
	public List<Triplet> generateForFactor(String factor, Set<String> mods, Set<String> tars,
		GeneProvider loader)
	{
		if (loader.get(factor) == null) return Collections.emptyList();

		// Don't use the modulators that are also targets
		mods.removeAll(tars);
		mods.remove(factor);
		tars.remove(factor);

		List<Triplet> trips = new ArrayList<>();

		mods.stream().filter(m -> loader.get(m) != null).forEach(m ->
			tars.stream().filter(t -> loader.get(t) != null).forEach(t ->
				trips.add(new Triplet(loader.get(m), loader.get(factor), loader.get(t)))));

		return trips.stream().filter(Triplet::hasAllGenes).collect(Collectors.toList());
	}

	/**
	 * Generates custom triplets for the given factor.
	 */
	public List<Triplet> generateForFactorAndModulator(String factor, String modulator, Set<String> tars,
		GeneProvider loader)
	{
		if (loader.get(factor) == null) return Collections.emptyList();

		// Don't use the modulators that are also targets
		tars.remove(factor);
		tars.remove(modulator);

		List<Triplet> trips = new ArrayList<>();

		tars.stream().filter(t -> loader.get(t) != null).forEach(t ->
				trips.add(new Triplet(loader.get(modulator), loader.get(factor), loader.get(t))));

		return trips.stream().filter(Triplet::hasAllGenes).collect(Collectors.toList());
	}
}
