/*
 * Copyright (c) 2014. Mount Sinai School of Medicine
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

package org.bdgenomics.avocado.algorithms.hmm

import scala.math.Ordering
import org.bdgenomics.adam.rich.RichADAMRecord

/**
 * Haplotype generated from HMM alignment.
 *
 * @param sequence String representing haplotype alignment.
 */
class Haplotype(val sequence: String, region: Seq[RichADAMRecord], hmm: HMMAligner = new HMMAligner, val reference: String = "") {

  lazy val referenceAlignment = hmm.alignSequences(reference, sequence, null)
  lazy val hasVariants = referenceAlignment.hasVariants

  lazy val perReadLikelihoods: Seq[Double] = region.map(read => {
    try {
      val alignment = HMMAligner.align(sequence, read.getSequence.toString, null)
      alignment.likelihood + alignment.prior
    } catch {
      case _: Throwable => {
        0.0
      }
    }
  })

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
    if (h1.readsLikelihood < h2.readsLikelihood) {
      -1
    } else {
      1
    }
  }
}
