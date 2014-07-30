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

import scala.math._
import scala.Ordering

object HaplotypePair {

  /**
   * Exponentiates two numbers by base of 10, adds together, takes base 10 log, and returns.
   *
   * @param x1 First digit to sum.
   * @param x2 Second digit to sum.
   * @return Sum of two digits after exponentation and logarithm.
   *
   * @see approxLogSumExp10
   */
  def exactLogSumExp10(x1: Double, x2: Double): Double = {
    log10(pow(10.0, x1) + pow(10.0, x2))
  }

  /**
   * Exponentiates two numbers by base of 10, adds together, takes base 10 log, and returns.
   *
   * @param x1 First digit to sum.
   * @param x2 Second digit to sum.
   * @return Sum of two digits after exponentation and logarithm.
   *
   * @see exactLogSumExp10
   */
  def approxLogSumExp10(x1: Double, x2: Double): Double = {
    exactLogSumExp10(x1, x2)
  }

}

/**
 * Class for a pairing of two haplotypes.
 *
 * @param haplotype1 First haplotype of pair.
 * @param haplotype2 Second haplotype of pair.
 */
class HaplotypePair(val haplotype1: Haplotype,
                    val haplotype2: Haplotype,
                    hmm: HMMAligner = new HMMAligner) {

  lazy val hasVariants = haplotype1.hasVariants || haplotype2.hasVariants
  lazy val (pairLikelihood, readAssignments) = scorePairLikelihood
  lazy val variantCount = haplotype1.variantCount + haplotype2.variantCount

  override def toString(): String = {
    haplotype1.sequence + ", " + haplotype1.hasVariants + "\n" +
      haplotype2.sequence + ", " + haplotype2.hasVariants + "\n" +
      ("%1.3f" format pairLikelihood) + ", " + hasVariants + "\n" +
      (0 until readAssignments.length).map(i => {
        readAssignments(i) + " <- " +
          haplotype1.perReadLikelihoods(i) + ", " +
          haplotype2.perReadLikelihoods(i) + ", Δ = " +
          (haplotype1.perReadLikelihoods(i) - haplotype2.perReadLikelihoods(i))
      }).reduce(_ + "\n" + _)
  }

  /**
   * Scores likelihood of two paired haplotypes and their alignment.
   *
   * @return Phred scaled likelihood.
   */
  def scorePairLikelihood: (Double, Array[Int]) = {
    val readLikelihoods = Array(haplotype1.perReadLikelihoods.toArray, haplotype2.perReadLikelihoods.toArray)
    val haplotypeLengths = Array(haplotype1.sequence.length.toLong, haplotype2.sequence.length.toLong)

    // run argmax
    // TODO: expose parameters for argmax fit
    HaplotypePairArgmax(haplotypeLengths, readLikelihoods, 10, 0.05, 0.5)
  }

  val sequences: Set[String] = Set(haplotype1.sequence, haplotype2.sequence)

  def sameSequences(other: HaplotypePair): Boolean = {
    sequences == other.sequences
  }

}

/**
 * Haplotype pairs are ordered by increasing pairwise likelihood, assuming
 *  they come from the same read group.
 */
object HaplotypePairOrdering extends Ordering[HaplotypePair] {

  /**
   * Compares two haplotype pairs. Returns 0 if the two pairs contain the same sequences. Else, sorts by pair
   * likelihood. If the likelihoods are the same, we break ties by looking at the count of variants in the pair
   * vs. the reference.
   *
   * @param pair1 First haplotype pair to compare.
   * @param pair2 Second haplotype pair to compare.
   * @return Comparison of haplotype pairs.
   */
  def compare(pair1: HaplotypePair, pair2: HaplotypePair): Int = {
    if (pair1.sameSequences(pair2)) {
      0
    } else if (pair1.pairLikelihood < pair2.pairLikelihood) {
      -1
    } else if (pair1.pairLikelihood > pair2.pairLikelihood) {
      1
    } else if (pair1.variantCount >= pair2.variantCount) {
      -1
    } else {
      1
    }
  }
}
