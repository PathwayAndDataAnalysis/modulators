package org.panda.gem.run;

import org.panda.gem.*;
import org.panda.gem.resource.CustomTripletMaker;
import org.panda.gem.resource.PCTripletMaker;
import org.panda.gem.resource.TCGAExpressionLoader;
import org.panda.utility.Kronometre;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.CorrelationSIFGenerator;
import org.panda.utility.statistics.Binomial;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * For a transcription factor, this class first generates its triplets using Pathway Commons, then tests them on 29 TCGA
 * RNAseq datasets, and writes resulting triplets in the results folder.
 *
 * @author Ozgun Babur
 */
public class TCGAPCRunner
{
//	static String[] codes = new String[]{"BRCA", "KIPAN", "GBMLGG", "LUAD", "THCA", "HNSC", "LUSC", "PRAD",
//		"SKCM", "COADREAD", "BLCA", "LIHC", "CESC", "OV", "UCEC", "STES", "ESCA", "PCPG", "PAAD", "LAML", "TGCT",
//		"THYM", "MESO", "UVM", "ACC", "UCS", "CHOL", "DLBC"};
	static String[] codes = new String[]{"BRCA"};

	static Map<String, Set<String>> selectedSubtypes = new HashMap<>();
	static
	{
//		selectedSubtypes.put("BRCA", new HashSet<>(Arrays.asList("LumA", "LumB", "Her2")));
	}

	static String outDir = "/home/ozgun/Analyses/GEM-runs/TCGA-PC/";

	public static void main(String[] args) throws IOException
	{
		Kronometre k = new Kronometre();
		String factor = "KLF4";
		String modulator = "ESR1";
		runInAllStudies(factor, modulator);
		integrate(factor + (modulator == null ? "" : "-" + modulator));
		k.print();
	}

	public static void runInAllStudies(String factor, String modulator) throws IOException
	{
		Map<String, Set<String>> subsets = readSubsets();

		Arrays.asList(codes).stream().forEach(code -> {
			try
			{
				System.out.println("code = " + code);
				TCGAExpressionLoader loader = new TCGAExpressionLoader("/home/ozgun/Data/TCGA/" + code,
					subsets.get(code));

				PCTripletMaker maker = new PCTripletMaker();
//				CustomTripletMaker maker = new CustomTripletMaker();
//				Set<String> targets = Files.lines(Paths.get("/home/ozgun/Documents/ESR1-responsive-genes.txt")).filter(l -> !l.isEmpty()).collect(Collectors.toSet());

				List<Triplet> trips = modulator == null ? maker.generateForFactor(factor, loader) :
					maker.generateForFactorAndModulator(factor, modulator, loader);

//				List<Triplet> trips = modulator == null ? maker.generateForFactor(factor, Collections.singleton(modulator), targets, loader) :
//					maker.generateForFactorAndModulator(factor, modulator, targets, loader);

				System.out.println("Triplet initial size = " + trips.size());

				trips = Selector.selectSignificantAndCategorized(trips, 0.1, 0.05);
				System.out.println("Triplet significant size = " + trips.size());
				write(trips, factor + (modulator == null ? "" : "-" + modulator), code);
			}
			catch(IOException e){throw new RuntimeException(e);}
		});
	}

	private static Map<String, Set<String>> readSubsets() throws IOException
	{
		if (selectedSubtypes != null && !selectedSubtypes.isEmpty())
		{
			Map<String, Set<String>> map = new HashMap<>();

			Files.lines(Paths.get("/media/ozgun/6TB/TCGA-pancan/pancan_samples.txt")).skip(1)
				.map(l -> l.split("\t")).filter(t -> selectedSubtypes.keySet().contains(t[2])).forEach(t ->
			{
				if (!map.containsKey(t[2])) map.put(t[2], new HashSet<>());
				if (selectedSubtypes.get(t[2]).contains(t[3])) map.get(t[2]).add(t[1]);
			});

			return map;
		}
		return Collections.emptyMap();
	}

	private static void write(List<Triplet> trips, String runName, String add) throws IOException
	{
		if (!trips.isEmpty())
		{
			Files.createDirectories(Paths.get(outDir + runName));

			String out = outDir + runName + File.separator + runName + "_" + add;
			Triplet.write(trips, out + ".txt");
			ModPrint mp = new ModPrint();
			mp.generateGEMPlot(trips, out + ".svg");
			if (!add.contains("recurrent")) writeModulatorCorrelations(trips, out);
		}
	}

	public static void integrate(String run) throws IOException
	{
		if (!Files.exists(Paths.get(outDir + run))) return;

		Map<Triplet, Integer> cnt = new HashMap<>();
		Files.newDirectoryStream(Paths.get(outDir + run)).forEach(p -> {
			if (p.toString().endsWith(".txt"))
			{
				try
				{
					Triplet.load(p.toString()).forEach(t -> cnt.put(t, cnt.containsKey(t) ? cnt.get(t) + 1 : 1));
				}
				catch (IOException e){throw new RuntimeException(e);}
			}
		});

		List<Triplet> trips;
		int rec = 2;

		do
		{
			int i = rec;
			trips = cnt.keySet().stream().filter(t -> cnt.get(t) >= i).collect(Collectors.toList());
			write(trips, run, "recurrent_" + rec);
			rec++;
		} while (!trips.isEmpty());


		write(cnt.keySet().stream().filter(t -> cnt.get(t) >= 3).collect(Collectors.toList()), run,
			"recurrent_" + 3);
	}

	public static void writeModulatorCorrelations(List<Triplet> trips, String filenameWithoutExtension)
	{
		Map<String, double[]> geneMap = new HashMap<>();
		trips.stream().map(t -> t.M).distinct().forEach(g -> geneMap.put(g.symbol, g.vals));
		CorrelationSIFGenerator csg = new CorrelationSIFGenerator(geneMap);
		csg.setMinCorrelation(0.3);

		Map<String, Map<ModulationCategory, Long>> cnt = trips.stream().collect(
			Collectors.groupingBy(t -> t.M.symbol, Collectors.groupingBy(t -> t.cat, Collectors.counting())));

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10},
			new Color[]{new Color(255, 200, 200), Color.WHITE, new Color(200, 255, 200)});

		cnt.keySet().forEach(g -> {
			int e = enhancerCnt(g, cnt);
			int a = attenuateCnt(g, cnt);
			double p = Binomial.getPval(a, e);
			double v = -Math.log(p) / Math.log(2);
			if (a > e) v = -v;

			csg.addNodeColor(g, vtc.getColor(v));
			csg.addNodeTooltip(g, String.valueOf(p));
		});

		csg.write(filenameWithoutExtension);
	}

	private static int enhancerCnt(String gene, Map<String, Map<ModulationCategory, Long>> cnt)
	{
		int c = 0;
		if (cnt.get(gene).containsKey(ModulationCategory.ENHANCES_ACTIVATION))
			c += cnt.get(gene).get(ModulationCategory.ENHANCES_ACTIVATION).intValue();
		if (cnt.get(gene).containsKey(ModulationCategory.ENHANCES_INHIBITION))
			c += cnt.get(gene).get(ModulationCategory.ENHANCES_INHIBITION).intValue();
		return c;
	}

	private static int attenuateCnt(String gene, Map<String, Map<ModulationCategory, Long>> cnt)
	{
		int c = 0;
		if (cnt.get(gene).containsKey(ModulationCategory.ATTENUATES_ACTIVATION))
			c += cnt.get(gene).get(ModulationCategory.ATTENUATES_ACTIVATION).intValue();
		if (cnt.get(gene).containsKey(ModulationCategory.ATTENUATES_INHIBITION))
			c += cnt.get(gene).get(ModulationCategory.ATTENUATES_INHIBITION).intValue();
		return c;
	}
}
