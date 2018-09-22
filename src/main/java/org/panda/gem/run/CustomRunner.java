package org.panda.gem.run;

import org.panda.gem.ModPrint;
import org.panda.gem.Selector;
import org.panda.gem.Triplet;
import org.panda.gem.resource.CustomTripletMaker;
import org.panda.gem.resource.GeneProvider;
import org.panda.gem.resource.SimpleFileExpressionLoader;
import org.panda.gem.resource.TCGAExpressionLoader;
import org.panda.utility.FileUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.panda.gem.run.SimpleRunner.TARGETS_FILE;

/**
 * This is for running GEM on custom set of modulators and targets, and expression dataset.
 *
 * Format of the custom expression file: Tab-delimited text file, first row is the header, first column is the gene
 * symbols, all other cells contain expression values. Each row (other than first) corresponds to a gene and each column
 * (other than first) corresponds to a sample in the dataset.
 *
 * Format of modulator and target candidate files: Gene symbols, one per line, no header.
 *
 * @author Ozgun Babur
 */
public class CustomRunner
{
	/**
	 * Name of the parameters file.
	 */
	private static final String PARAMETER_FILENAME = "parameters.txt";

	/**
	 * The directory that has the parameters.txt file.
	 */
	private String inputDirectory;

	/**
	 * Name of the custom expression file.
	 */
	private String customExpressionFile;

	/**
	 * The directory that contains TCGA datasets that are downloaded from Broad Firehose using the BroadDownloader in
	 * the "resource" project.
	 */
	private String tcgaDirectory;

	/**
	 * Study code of the TCGA disease.
	 */
	private String tcgaStudy;

	/**
	 * Set of the TCGA subtypes desired to be used by the analysis.
	 */
	private Set<String> tcgaSubtypes;

	/**
	 * Set of the modulators.
	 */
	private Set<String> modulators;

	/**
	 * Set of the targets.
	 */
	private Set<String> targets;

	/**
	 * The transcription factor to use.
	 */
	private String factor;

	/**
	 * The false discovery rate that will be applied on selecting significant gamma.
	 */
	private double fdrThr = 0.1;

	/**
	 * The p-value threshold to use while determining the category of a significant triplet.
	 */
	private double categoryPvalThr = 0.05;

	/**
	 * Name of the result triplet text file.
	 */
	private String tripletFilename = "results.txt";

	/**
	 * Name of the result SVG file that shows distribution of triplet categories for each modulator separately.
	 */
	private String svgFilename = "results.svg";


	public CustomRunner(String inputDirectory) throws IOException
	{
		this.inputDirectory = inputDirectory;
		this.modulators = new HashSet<>();
		this.targets = new HashSet<>();
		this.tcgaSubtypes = new HashSet<>();
		this.tripletFilename = inputDirectory + File.separator + this.tripletFilename;
		this.svgFilename = inputDirectory + File.separator + this.svgFilename;

		readParameters();
	}

	private void readParameters() throws IOException
	{
		Files.lines(Paths.get(inputDirectory + File.separator + PARAMETER_FILENAME)).filter(l -> !l.trim().startsWith("#"))
			.filter(l -> !l.isEmpty())
			.map(l -> new String[]{l.substring(0, l.indexOf("=")).trim(), l.substring(l.indexOf("=") + 1).trim()})
			.forEach(t -> setParameter(t[0], t[1]));
	}

	private void setParameter(String key, String value)
	{
		Parameter param = Parameter.findEnum(key);

		if (param == null) throw new RuntimeException("Unknown parameter: " + key);

		param.reader.read(value, this);
	}

	public void addModulator(String modulator)
	{
		this.modulators.add(modulator);
	}

	public void addTarget(String target)
	{
		this.targets.add(target);
	}

	protected String getFilename(String param)
	{
		if (param.startsWith(File.separator)) return param;
		return inputDirectory + File.separator + param;
	}

	interface ParameterReader
	{
		void read(String value, CustomRunner cr);
	}

	enum Parameter
	{
		MODULATOR((value, cr) -> cr.addModulator(value), "Modulator",
			"Gene symbol of a modulator to use in the analysis"),
		TARGET((value, cr) -> cr.addTarget(value), "Target", "Gene symbol of a target to use in the analysis."),
		FACTOR((value, cr) -> cr.factor = value, "Factor", "The gene symbol of the transcription factor."),
		MODULATORS_FILE((value, cr) -> cr.modulators.addAll(
			FileUtil.getLinesInSet(cr.getFilename(value))),
			"Modulators filename", "The file has to contain a gene symbol per line."),
		TARGETS_FILE((value, cr) -> cr.targets.addAll(FileUtil.getLinesInSet(cr.getFilename(value))),
			"Targets filename", "The file has to contain a gene symbol per line."),
		CUSTOM_EXPRESSION_FILE((value, cr) -> cr.customExpressionFile = cr.getFilename(value),
			"Name of the custom expression file",
			"First column has gene symbols, first row has sample names, values are in columns."),
		TCGA_DIRECTORY((value, cr) -> cr.tcgaDirectory = cr.getFilename(value), "TCGA data directory",
			"The directory where TCGA data is downloaded using the BroadDownloader class in the \"resource\" project."),
		TCGA_STUDY((value, cr) -> cr.tcgaStudy = value, "TCGA study code", "The disease code toget expression from."),
		USE_SUBTYPE((value, cr) -> cr.tcgaSubtypes.add(value), "Subtype",
			"Select one or more pre-defined subtypes to limit the analysis."),
		FDR_THR ((value, cr) -> cr.fdrThr = Double.valueOf(value), "FDR threshold",
			"The false discovery cutoff to use for determining significant modulations."),
		CATEGORY_PVAL_THR((value, cr) -> cr.categoryPvalThr = Double.valueOf(value), "Category p-value threhsold",
			"For using while determining category of a significant triplet. Note that some significant triplets may " +
				"fail to get a category. In that case they are excluded from results."),
		RESULT_TRIPLET_FILENAME((value, cr) -> cr.tripletFilename = cr.getFilename(value), "Result triplet filename",
			"Overrides the default name"),
		RESULT_SVG_FILENAME((value, cr) -> cr.svgFilename= cr.getFilename(value), "Result SVG filename",
			"Overrides the default name"),
		;

		ParameterReader reader;
		String title;
		String info;

		Parameter(ParameterReader reader, String title, String info)
		{
			this.reader = reader;
			this.title = title;
			this.info = info;
		}

		String getText()
		{
			return toString().toLowerCase().replaceAll("_", "-");
		}

		static Parameter findEnum(String text)
		{
			if (text == null) return null;

			text = text.trim();
			for (Parameter parameter : values())
			{
				if (parameter.getText().equals(text)) return parameter;
			}
			return null;
		}
	}

	public void run() throws IOException
	{
		GeneProvider loader = null;

		if (customExpressionFile != null)
		{
			// Load expression data
			loader = new SimpleFileExpressionLoader(customExpressionFile);
		}
		else if (tcgaDirectory != null && tcgaStudy != null)
		{
			Set<String> subsets = readSubsets();

			loader = new TCGAExpressionLoader(tcgaDirectory + File.separator + tcgaStudy, subsets);
		}

		// Prepare triplets using the custom modulators and targets sets.
		CustomTripletMaker maker = new CustomTripletMaker();
		List<Triplet> trips = maker.generateForFactor(factor, modulators, targets, loader);
		System.out.println("Size of triplets tested      = " + trips.size());

		// Select significant triplets and determine modulation categories
		trips = Selector.selectSignificantAndCategorized(trips, fdrThr, categoryPvalThr);
		System.out.println("Size of significant triplets = " + trips.size());

		// Write result triplets
		Triplet.write(trips, tripletFilename);

		// Draw the result graphic
		ModPrint mp = new ModPrint();
		mp.generateGEMPlot(trips, svgFilename);
	}

	private Set<String> readSubsets() throws IOException
	{
		if (!tcgaSubtypes.isEmpty())
		{
			Set<String> samples = new HashSet<>();

			Files.lines(Paths.get(tcgaDirectory + "/pancan_samples.txt")).skip(1)
				.map(l -> l.split("\t")).filter(t -> t[2].equals(tcgaStudy)).forEach(t ->
			{
				if (tcgaSubtypes.contains(t[3])) samples.add(t[1]);
			});

			return samples;
		}
		return null;
	}


	/**
	 * Example run of GEM
	 */
	public static void main(String[] args) throws IOException
	{
		List<String> dirs = FileUtil.getSubdirectoriesContaining(args[0], "parameters.txt");

		for (String dir : dirs)
		{
			System.out.println("directory = " + dir);
			CustomRunner cr = new CustomRunner(dir);
			cr.run();
		}
	}
}
