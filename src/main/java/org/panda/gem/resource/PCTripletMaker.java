package org.panda.gem.resource;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.gem.Triplet;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.GraphList;
import org.panda.utility.graph.UndirectedGraph;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Generates triplets (hypotheses) using Pathway Commons SIF graph. Modulator candidates are the neighbors in the
 * in-complex-with and interacts-with graphs, and upstream neighbors in the controls-state-change-of graphs.
 *
 * @author Ozgun Babur
 */
public class PCTripletMaker
{
	/**
	 * F -> T graph
	 */
	DirectedGraph expGraph;

	/**
	 * M - F PPI graph
	 */
	GraphList ppiGraph;

	/**
	 * M -> F state-change graph
	 */
	DirectedGraph stcGraph;

	public PCTripletMaker()
	{
		PathwayCommons pc = new PathwayCommons();
		expGraph = (DirectedGraph) pc.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF);
		ppiGraph = (GraphList) pc.getGraph(SIFEnum.IN_COMPLEX_WITH, SIFEnum.INTERACTS_WITH);
		stcGraph = (DirectedGraph) pc.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
	}

	/**
	 * Generates triplets from Pathway Commons for the given factor.
	 */
	public List<Triplet> generateForFactor(String factor, GeneProvider loader)
	{
		if (loader.get(factor) == null) return Collections.emptyList();

		Set<String> mods = new HashSet<>(ppiGraph.getNeighbors(factor));
		mods.addAll(stcGraph.getUpstream(factor));
		Set<String> tars = new HashSet<>(expGraph.getDownstream(factor));
		if (factor.equals("MYC")) tars.addAll(MYC_TARGETS);

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
	 * Generates triplets from Pathway Commons for the given factor.
	 */
	public List<Triplet> generateForFactorAndModulator(String factor, String modulator, GeneProvider loader)
	{
		if (loader.get(factor) == null) return Collections.emptyList();

		Set<String> tars = new HashSet<>(expGraph.getDownstream(factor));
		if (factor.equals("MYC")) tars.addAll(MYC_TARGETS);

		// Don't use the modulators that are also targets
		tars.remove(factor);

		List<Triplet> trips = new ArrayList<>();

		tars.stream().filter(t -> loader.get(t) != null).forEach(t ->
				trips.add(new Triplet(loader.get(modulator), loader.get(factor), loader.get(t))));

		return trips.stream().filter(Triplet::hasAllGenes).collect(Collectors.toList());
	}

	private static Set<String> MYC_TARGETS = new HashSet<>(Arrays.asList((
		"E2F3\n" +
		"PFKM\n" +
		"PAPPA-AS1\n" +
		"EP300\n" +
		"CLN5\n" +
		"POLR3D\n" +
		"TP53\n" +
		"SERPINI1\n" +
		"SUPT7L\n" +
		"SUPT3H\n" +
		"MYCT1\n" +
		"MYC\n" +
		"FOSL1\n" +
		"CCNB1\n" +
		"CCND2\n" +
		"MTDH\n" +
		"LDHA\n" +
		"PRDX3\n" +
		"LIN28B\n" +
		"MMP9\n" +
		"IREB2\n" +
		"IRG1\n" +
		"ODC1\n" +
		"PTMAP4\n" +
		"PIM1\n" +
		"KAT2A\n" +
		"KAT5\n" +
		"GPAM\n" +
		"TFRC\n" +
		"PEG10\n" +
		"DDX18\n" +
		"RUVBL2\n" +
		"RUVBL1\n" +
		"TERT\n" +
		"ID2\n" +
		"SMAD4\n" +
		"HSPD1\n" +
		"HSPA4\n" +
		"PTMAP7\n" +
		"PTMA\n" +
		"HSP90AA1\n" +
		"SHMT1\n" +
		"SMAD3\n" +
		"PMAIP1\n" +
		"DFFB\n" +
		"NBN\n" +
		"MTA1\n" +
		"MINA\n" +
		"GAPDH\n" +
		"TAF10\n" +
		"TAF9\n" +
		"TAF12\n" +
		"GAPDHP44\n" +
		"UBTF\n" +
		"EIF2S1\n" +
		"EIF4E\n" +
		"EIF4A1\n" +
		"HMGA1\n" +
		"TAF4B\n" +
		"ENO1\n" +
		"SLC2A1\n" +
		"SLC25A21\n" +
		"HUWE1\n" +
		"TRRAP\n" +
		"TK1\n" +
		"CDK4\n" +
		"NDUFAF2\n" +
		"NPM1\n" +
		"CDC25A\n" +
		"CDCA7\n" +
		"EIF4G1\n" +
		"RMRP\n" +
		"RCC1\n" +
		"COMMD3-BMI1\n" +
		"MAX\n" +
		"BCAT1\n" +
		"RPL11\n" +
		"SNAI1\n" +
		"ACTL6A\n" +
		"LONP1\n" +
		"BMI1\n" +
		"DNAJC5\n" +
		"CALD1\n" +
		"CAD\n" +
		"BAX\n" +
		"NME1-NME2\n" +
		"NME2\n" +
		"NME1\n" +
		"KIR3DL1\n" +
		"NCL\n" +
		"PDCD10\n" +
		"CREBBP\n" +
		"BIRC5\n" +
		"ARTN").split("\n")));
}
