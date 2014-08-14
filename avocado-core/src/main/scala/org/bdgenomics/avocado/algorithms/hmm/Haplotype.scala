/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.avocado.algorithms.hmm

import scala.math.Ordering
import org.bdgenomics.adam.rich.RichAlignmentRecord
import org.bdgenomics.adam.util.PhredUtils
import scala.math._

/**
 * Haplotype generated from HMM alignment.
 *
 * @param sequence String representing haplotype alignment.
 */
class Haplotype(val sequence: String,
                region: Seq[RichAlignmentRecord],
                hmm: HMMAligner = new HMMAligner,
                val reference: String = "",
                readAlignerConfig: TransitionMatrixConfiguration = TransitionMatrixConfiguration()) {

  def reg: String = {
    val start = region.map(_.getStart).min
    val end = region.map(_.getEnd).max
    val name = region.head.getContig.getContigName
    name + ", " + start + ", " + end
  }

  val fastAligner = new FastAligner(sequence, region.head.record.getSequence.length)
  val exactAligner = new ExactReadAligner(sequence, region.head.record.getSequence.length)

  assert(reference.length > 0, "Reference has length 0 on " + reg + ".")
  assert(sequence.length > 0, "Haplotype has length 0 on " + reg + ".")
  lazy val referenceAlignment = hmm.alignSequences(reference, sequence, null)
  assert(referenceAlignment.hasVariants == (sequence != reference), "HMM calls variant, but sequence matches reference. " + referenceAlignment)
  lazy val hasVariants = referenceAlignment.hasVariants
  lazy val variantCount = referenceAlignment.alignmentStateSequence
    .filter(s => s == 'X' || s == 'I' || s == 'D').length

  lazy val perReadAlignments: Seq[Alignment] = {
    region.map(read => {
      try {
        fastAligner.alignSequences(sequence, read.getSequence.toString, read.getQual.toString)
      } catch {
        case _: Throwable => {
          exactAligner.alignSequences(sequence, read.getSequence.toString, read.getQual.toString)
        }
      }
    })
  }
  lazy val perReadLikelihoods: Seq[Double] = perReadAlignments.map(_.likelihood)
  lazy val readsLikelihood = perReadLikelihoods.sum

  override def toString(): String = {
    sequence
  }

}

/**
 * Haplotypes are ordered by increasing reads likelihood, assuming they
 * come from the same group of reads.
 */
object HaplotypeOrdering extends Ordering[Haplotype] {

  /**
   * Compares two haplotypes. Returns (-1, 0, 1) if h1 has (lower, same, higher) read
   * likelihood than h2. 0 is only returned if they have the same sequence AND likelihood
   *
   * @param h1 First haplotype to compare.
   * @param h2 Second haplotype to compare.
   * @return Ordering info for haplotypes.
   */
  def compare(h1: Haplotype, h2: Haplotype): Int = {
    if (h1.sequence == h2.sequence) {
      h1.readsLikelihood.compare(h2.readsLikelihood)
    } else if (h1.readsLikelihood < h2.readsLikelihood) {
      -1
    } else {
      1
    }
  }
}
