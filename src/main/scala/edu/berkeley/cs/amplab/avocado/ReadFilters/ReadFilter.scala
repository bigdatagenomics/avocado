package edu.berkeley.cs.amplab.avocado.filters.reads

import spark.{RDD,SparkContext}
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord

/**
 * Trait for filtering reads.
 */
trait ReadFilter {

  /**
   * Method signature for filter operation.
   *
   * @param[in] reads An RDD containing reference oriented stacks of nucleotides.
   * @return An RDD containing lists of pileups.
   */
  def filter (reads: RDD [(void, ADAMRecord)]): RDD [(Any, List[ADAMRecord])]
}
