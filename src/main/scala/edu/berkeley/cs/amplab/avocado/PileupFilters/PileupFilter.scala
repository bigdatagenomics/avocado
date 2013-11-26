/*
 * Copyright (c) 2013. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.berkeley.cs.amplab.avocado.filters.pileup

import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.ADAMPileup
import edu.berkeley.cs.amplab.avocado.Avocado

/**
 * Abstract class for filtering pileups. 
 */
abstract class PileupFilter extends Serializable with Logging {

  val filterName: String
  val coverage = Avocado.coverage
  val reducePartitions = Avocado.reducePartitions

  /**
   * Method signature for filter operation.
   *
   * @param[in] pileups An RDD containing reference oriented stacks of nucleotides.
   * @return An RDD containing lists of pileups.
   */
  def filter (pileups: RDD [ADAMPileup]): RDD [ADAMPileup]
}
