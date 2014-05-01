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
class HaplotypePair(val haplotype1: Haplotype, val haplotype2: Haplotype, hmm: HMMAligner = new HMMAligner) {

  lazy val hasVariants = haplotype1.hasVariants || haplotype2.hasVariants
  lazy val pairLikelihood = scorePairLikelihood

  override def toString(): String = {
    haplotype1.sequence + ", " + haplotype2.sequence + ", " + ("%1.3f" format pairLikelihood)
  }

  /**
   * Scores likelihood of two paired haplotypes and their alignment.
   *
   * @return Phred scaled likelihood.
   */
  def scorePairLikelihood: Double = {
    val readsProb = haplotype1.perReadLikelihoods.zip(haplotype1.perReadLikelihoods).map(scores => HaplotypePair.exactLogSumExp10(scores._1, scores._2) - log10(2.0)).sum
    val alignment = hmm.alignSequences(haplotype2.sequence, haplotype1.sequence, null)
    readsProb + alignment.prior
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
   * Compares two haplotype pairs. Returns (-1, 0, 1) if first pair has (lower, same, higher)
   * pairwise likelihood.
   *
   * @param pair1 First haplotype pair to compare.
   * @param pair2 Second haplotype pair to compare.
   * @return Comparison of haplotype pairs.
   */
  def compare(pair1: HaplotypePair, pair2: HaplotypePair): Int = {
    if (pair1.pairLikelihood < pair2.pairLikelihood) {
      -1
    } else {
      1
    }
  }
}