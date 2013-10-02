package edu.berkeley.cs.amplab.avocado.calls.pileup

import spark.{RDD,SparkContext}
import edu.berkeley.cs.amplab.adam.util.{Pileup,PileupTraversable}
import edu.berkeley.cs.amplab.adam.avro.ADAMVariant

/**
 * Trait for filtering pileups. 
 */
trait PileupCall {

  /**
   * Method signature for variant calling operation.
   *
   * @param[in] pileupGroups An RDD containing lists of pileups.
   * @return An RDD containing called variants.
   */
  def call (pileupGroups: RDD [(void, List[Pileup])]): RDD [ADAMVariant])] 
}

