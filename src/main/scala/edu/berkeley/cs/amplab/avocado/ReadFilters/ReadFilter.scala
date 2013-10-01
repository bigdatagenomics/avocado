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
   * @param[in] reads An RDD containing reads.
   * @return An RDD containing lists of reads.
   */
  def filter (reads: RDD [(void, ADAMRecord)]): RDD [(Any, List[ADAMRecord])]
}
