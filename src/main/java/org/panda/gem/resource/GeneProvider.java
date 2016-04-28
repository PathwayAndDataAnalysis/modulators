package org.panda.gem.resource;

import org.panda.gem.Gene;

/**
 * Interface for any resource that can provide a gene to use in GEM analysis.
 * @author Ozgun Babur
 */
public interface GeneProvider
{
	/**
	 * Should return (create if not already created) the gene with its expression array initialized.
	 * @param symbol HGNC symbol of the gene
	 */
	Gene get(String symbol);
}
