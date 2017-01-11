package org.panda.gem.run;

import org.panda.gem.ModPrint;
import org.panda.gem.Selector;
import org.panda.gem.Triplet;
import org.panda.gem.resource.CustomTripletMaker;
import org.panda.gem.resource.SimpleFileExpressionLoader;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A sample execution of GEM using a custom expression data file, for a transcription factor, with custom set of
 * modulator candidates, and custom target candidates.
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
	static final String DIR = "/home/babur/Documents/expO/";

	/**
	 * Name of the expression file.
	 */
	static final String EXPRESSION_FILE = DIR + "GSE2109_series_matrix.txt";

	/**
	 * Name of the modulators file.
	 */
	static final String MODULATORS_FILE = DIR + "modulators.txt";

	/**
	 * Name of the targets file.
	 */
	static final String TARGETS_FILE = DIR + "targets.txt";

	/**
	 * Gene symbol fo the transcription factor in focus.
	 */
	static final String TF_SYMBOL = "AR";

	/**
	 * The false discovery rate that will be applied on selecting significant gamma.
	 */
	static final double FDR_THR = 0.05;

	/**
	 * The p-value threshold to use while determining the category of a significant triplet.
	 */
	static final double CATEG_PVAL_THR = 0.05;

	/**
	 * Name of the result triplet text file.
	 */
	static final String RESULT_TRIPLET_FILE = DIR + "GEM-result-" + TF_SYMBOL + ".txt";

	/**
	 * Name of the result SVG file that shows distribution of triplet categories for each modulator separately.
	 */
	static final String RESULT_SVG_FILE = DIR + "GEM-result-" + TF_SYMBOL + ".svg";

	/**
	 * Example run of GEM
	 */
	public static void main(String[] args) throws IOException
	{
		// Load expression data
		SimpleFileExpressionLoader loader = new SimpleFileExpressionLoader(EXPRESSION_FILE);

		// Prepare triplets using the custom modulators and targets sets.
		CustomTripletMaker maker = new CustomTripletMaker();
		Set<String> modulators = Files.lines(Paths.get(MODULATORS_FILE)).collect(Collectors.toSet());
		Set<String> targets = Files.lines(Paths.get(TARGETS_FILE)).collect(Collectors.toSet());
		List<Triplet> trips = maker.generateForFactor(TF_SYMBOL, modulators, targets, loader);

		// Select significant triplets and determine modulation categories
		trips = Selector.selectSignificantAndCategorized(trips, FDR_THR, CATEG_PVAL_THR);
		System.out.println("Size of significant triplets = " + trips.size());

		// Write result triplets
		Triplet.write(trips, RESULT_TRIPLET_FILE);

		// Draw the result graphic
		ModPrint mp = new ModPrint();
		mp.generateGEMPlot(trips, RESULT_SVG_FILE);
	}
}
