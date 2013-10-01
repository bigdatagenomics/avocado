package edu.berkeley.cs.amplab.avocado.filters.pileup

import spark.{RDD,SparkContext}
import edu.berkeley.cs.amplab.adam.util.{Pileup,PileupTraversable}
import edu.berkeley.cs.amplab.avocado.filters.pileup.PileupFilter

/**
 * Class implementing a filter that only returns Pileups that contain
 * a mismatch on any bases in the pileup.
 */
class PileupFilterOnMismatch extends PileupFilter {

  /**
   * Filter to only return pileups with mismatches.
   *
   * @param[in] pileups An RDD containing reference oriented stacks of nucleotides.
   * @return An RDD containing only pileups that contain at least one mismatch.
   */
  override def filter (pileups: RDD [(void, Pileup)]): RDD [(void, List[Pileup])] = {
    return pileups.filter (_.mismatches != List.empty)
                  .map (_ => (void, List (_)))
  }
}
