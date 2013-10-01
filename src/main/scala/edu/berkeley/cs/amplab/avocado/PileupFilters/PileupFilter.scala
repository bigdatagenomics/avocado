package edu.berkeley.cs.amplab.avocado.filters.pileup

import spark.{RDD,SparkContext}
import edu.berkeley.cs.amplab.adam.util.{Pileup,PileupTraversable}

/**
 * Trait for filtering pileups. 
 */
trait PileupFilter {

  /**
   * Method signature for filter operation.
   *
   * @param[in] pileups An RDD containing reference oriented stacks of nucleotides.
   * @return An RDD containing lists of pileups.
   */
  def filter (pileups: RDD [(void, Pileup)]): RDD [(void, List[Pileup])] 
}
