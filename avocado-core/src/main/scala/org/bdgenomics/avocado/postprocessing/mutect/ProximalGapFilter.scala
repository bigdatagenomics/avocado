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

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.adam.rdd.RegionRDDFunctions._
import org.bdgenomics.formats.avro.AlignmentRecord

import ClassifiedContext._

object ProximalGapFilter {

  /**
   * Different descriptions of Mutect have different descriptions of this filter --
   * it's not clear whether it's
   *   (>= 3 reads with either an insertion or a deletion)
   * or
   *   (>= 3 reads with an insertion) or (>= 3 reads with a deletion)
   *
   * The paper seems to suggest the latter, but other sources (includes slides I've found)
   * suggest the former.
   *
   * So I'm abstracting it out, with a simple implementation here that can be parameterized
   * later.
   */
  def hasIndel(rec: AlignmentRecord): Boolean =
    rec.getCigar.contains("I") || rec.getCigar.contains("D")
}

/**
 * ProximalGapFilter implements the HC filter described by (in the paper) the following
 * text:
 *   Remove false positives caused by nearby misaligned small insertion and deletion events.
 *   Reject candidate site if there are ≥ 3 reads with insertions in an 11-base-pair window
 *   centered on the candidate mutation or if there are ≥ 3 reads with deletions in the same
 *   11-base-pair window.
 *
 * See the notes to the ProximalGapFilter.hasIndel method, above, for more details on some
 * ambiguity in this statement.
 *
 * @param indelReadThreshold The number of reads passing the filter below which the variant is
 *                           retained (by default == 3)
 * @param windowDistance The distance around the variant to look for overlapping reads
 *                       (by default == 11)
 * @param gapFilter The function that determines whether or not a read has an indel.
 */
class ProximalGapFilter(val indelReadThreshold: Int = 3,
                        val windowDistance: Int = 11,
                        val gapFilter: AlignmentRecord => Boolean = ProximalGapFilter.hasIndel)
    extends MutectPostprocessor with Serializable {

  override def filter(variants: RDD[VariantContext],
                      tumorReads: RDD[Classified[AlignmentRecord]],
                      normalReads: RDD[Classified[AlignmentRecord]]): RDD[VariantContext] =
    variants
      .groupWithinRange(
        tumorReads.filterByClasses("retained").values(),
        windowDistance.toLong)
      .filter(_._2.count(gapFilter) < indelReadThreshold)
      .map(_._1)
}
