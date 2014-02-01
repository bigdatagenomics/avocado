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

package edu.berkeley.cs.amplab.avocado.calls.pileup

import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._
import edu.berkeley.cs.amplab.adam.rdd.AdamRDDFunctions
import edu.berkeley.cs.amplab.adam.models.{ADAMRod, ADAMVariantContext}
import edu.berkeley.cs.amplab.avocado.calls.VariantCall
import edu.berkeley.cs.amplab.avocado.Avocado
import org.apache.spark.rdd.RDD
import org.apache.spark.{SparkContext, Logging}

/**
 * Abstract class for calling variants on reads. 
 */
abstract class PileupCall extends VariantCall {

  /**
   * Method signature for variant calling operation.
   *
   * @param pileups An RDD of pileups.
   * @return An RDD containing called variants.
   */
  def callRods (pileups: RDD[ADAMRod]): RDD[ADAMVariantContext]

  /**
   * Converts reads to rods, then calls rod based variant caller.
   *
   * @param reads An RDD of reads.
   * @return An RDD containing called variants.
   */
  def call (reads: RDD[ADAMRecord]): RDD[ADAMVariantContext] = {
    val rods = reads.adamRecords2Rods()
    callRods(rods)
  }
  
  override def isReadCall (): Boolean = false
  override def isPileupCall (): Boolean = true
}
