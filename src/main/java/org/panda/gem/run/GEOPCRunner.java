package org.panda.gem.run;

import org.panda.gem.ModPrint;
import org.panda.gem.Selector;
import org.panda.gem.Triplet;
import org.panda.gem.resource.GEOFromFileLoader;
import org.panda.gem.resource.PCTripletMaker;

import java.io.IOException;
import java.util.List;

/**
 * A sample execution of GEM using the expO dataset in Gene Expression Omnibus database.
 *
 * @author Ozgun Babur
 */
public class GEOPCRunner
{
	static final String DIR = "/home/babur/Documents/expO/";

	/**
	 * This is the platform file for Affymetrix Human Genome U133 Plus 2.0 Array. Change it accordingly.
	 */
	static final String PLATFORM_FILE = DIR + "GPL570.txt";

	/**
	 * Name of the series file that contains expression values. This one is for expO dataset.
	 */
	static final String SERIES_FILE = DIR + "GSE2109_series_matrix.txt";

	/**
	 * Gene symbol fo the transcription factor in focus.
	 */
	static final String TF_SYMBOL = "CLK1";

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
		// Load GEO data
		GEOFromFileLoader loader = new GEOFromFileLoader(PLATFORM_FILE, SERIES_FILE);

		// Prepare triplets using Pathway Commons 2 version 8. In this case, most of the transcriptional targets come
		// from TRANSFAC, and most of the binding proteins come from IntAct and HPRD, while there are also many other
		// databases contributing.
		PCTripletMaker maker = new PCTripletMaker();
		List<Triplet> trips = maker.generateForFactor(TF_SYMBOL, loader);

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
