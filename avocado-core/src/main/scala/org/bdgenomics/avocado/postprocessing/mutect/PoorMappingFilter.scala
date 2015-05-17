/*
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.avocado.postprocessing.mutect

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.adam.rich.RichAlignmentRecord._
import org.bdgenomics.formats.avro.AlignmentRecord

import scala.math._

/**
 * Implements the Poor Mapping Filter, described in the paper with the following text:
 *
 *   Remove false positives caused by sequence similarity in the genome, leading to
 *   misplacement of reads. Two tests are used to identify such sites:
 *     (i) candidates are rejected if ≥ 50% of the reads in the tumor and normal
 *         samples have a mapping quality of zero (although reads with a mapping
 *         quality of zero are discarded in the short-read pre-processing (Supplementary
 *         Methods), this filter reconsiders those discarded reads); and
 *     (ii) candidates are rejected if they do not have at least a single observation
 *          of the mutant allele with a confident mapping (that is, mapping quality
 *          score ≥ 20).
 */
class PoorMappingFilter(val mapq0Fraction: Double = 0.5,
                        val confidentMappingQuality: Int = 20) extends MutectPostprocessor with Serializable {

  // TODO(twd): need to review this, not sure a mismatch is the right standard for "observes mutant allele."
  def observesMutantAllele(vc: VariantContext)(read: AlignmentRecord): Boolean = {
    read.isMismatchAtReferencePosition(vc.position).getOrElse(false)
  }

  override def filter(variants: RDD[VariantContext],
                      tumorReads: RDD[Classified[AlignmentRecord]],
                      normalReads: RDD[Classified[AlignmentRecord]]): RDD[VariantContext] = {

    val tumorMapqFractions: RDD[(VariantContext, (Int, Int))] = ???
    //      variants
    //      .groupByOverlap(tumorReads.values())
    //      .map(p => (p._1, (p._2.count(r => r.getMapq == 0), p._2.size)))

    val normalMapqFractions: RDD[(VariantContext, (Int, Int))] = ???
    //      variants
    //      .groupByOverlap(normalReads.values())
    //      .map(p => (p._1, (p._2.count(r => r.getMapq == 0), p._2.size)))

    // This captures the set of variants which _pass_ condition (i)
    val filter1 = tumorMapqFractions.join(normalMapqFractions)
      .filter {
        case (vc: VariantContext, counts: ((Int, Int), (Int, Int))) => {

          val zeroed = counts._1._1 + counts._2._1
          val total = counts._1._2 + counts._2._2
          val fraction = zeroed.toDouble / max(1, total)

          fraction < mapq0Fraction
        }
      }.map(_._1)

    // This captures the set of variants which _pass_ condition (ii)
    val filter2 = ???
    //      variants.groupByOverlap(tumorReads.filterByClasses("retained").values())
    //      .filter {
    //        case (vc: VariantContext, reads: Iterable[AlignmentRecord]) =>
    //          reads.exists(observesMutantAllele(vc))
    //      }.map(_._1)

    // And we return the intersection of the variants (i.e. that pass _both_ filters)
    filter1.intersection(filter2)
  }
}

