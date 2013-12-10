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

package edu.berkeley.cs.amplab.avocado.filters.reads

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import scala.math.{min,max}
import edu.berkeley.cs.amplab.avocado.calls.reads.ReadCallAssemblyPhaser
import edu.berkeley.cs.amplab.avocado.calls.pileup.PileupCallUnspecified
import edu.berkeley.cs.amplab.avocado.calls.VariantCall

/**
 * Enumeration for levels of complexity (high, medium, low).
 */
object MapComplexity extends Enumeration {
  val High, Low = Value
}

/**
 * Class that implements a complexity based read filter.
 */ 
class ReadFilterOnComplexity (coverageThresholdHigh: Int,
                              coverageThresholdLow: Int,
                              mappingQualityThreshold: Int) extends ReadFilter {

  val filterName = "ComplexityRegions"

  val stripe = 1000
  val overlap = 100

  /**
   * Calculates the complexity of a region. Uses the following heuristics:
   *  * Mapping Quality: What is the average mapping quality of this region?
   *  * Coverage: What is the average coverage over this region?
   * These values are compared against thresholds set in configuration.
   *
   * @param[in] reads A list of reads.
   * @return Tuple of region complexity and list of reads.
   */
  def scoreComplexity (segment: (Long, Seq[ADAMRecord])): (MapComplexity.Value, Seq[ADAMRecord]) = {
    val (segmentNumber, reads) = segment
    val segmentStart = segmentNumber * stripe
    val segmentEnd = (segmentNumber + 1) * stripe
    val complexity = reads.map (_.getMapq).reduce (_ + _) / reads.length
    val coverage = reads.map (v => min (v.getStart + v.getSequence.length, segmentEnd) - max (v.getStart, segmentStart))
                        .reduce (_ + _) / reads.length

    if (complexity < mappingQualityThreshold || (coverage > coverageThresholdHigh) || (coverage < coverageThresholdLow)) {
      return (MapComplexity.High, reads)
    } else {
      return (MapComplexity.Low, reads)
    }
  }

  /**
   * Simple complexity filter. Looks at average mapping quality and pileup depth.
   *
   * @param[in] reads An RDD containing reads.
   * @return An RDD containing lists of reads.
   */
  def filter (reads: RDD [ADAMRecord]): Map[VariantCall, RDD[ADAMRecord]] = {

    var callMap = Map[VariantCall, RDD[ADAMRecord]]()
    
    /* split into windows of reads
     * FIXME: add overlap on windows
     */
    val segments = reads.groupBy (v => v.getStart / stripe)
    
    // compute complexity value, group by complexity, and return
    val segmentsByComplexity = segments.map (scoreComplexity)

    // call high compexity regions with a read call
    val readCallRegions = segmentsByComplexity.filter (_._1 == MapComplexity.High)
      .flatMap (_._2)
    callMap += (new ReadCallAssemblyPhaser -> readCallRegions)

    // call low complexity regions with a pileup call
    val pileupCallRegions = segmentsByComplexity.filter (_._1 == MapComplexity.Low)
      .flatMap (_._2)
    callMap += (new PileupCallUnspecified -> pileupCallRegions)

    // return call map
    callMap
  }
}
