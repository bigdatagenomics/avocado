package org.bdgenomics.avocado.algorithms.hmm

import scala.collection.mutable.ArrayBuffer
import org.bdgenomics.adam.avro.ADAMRecord
import scala.math._
import scala.math.Ordering
import org.bdgenomics.adam.rich.RichADAMRecord

/**
 * Haplotype generated from HMM alignment.
 *
 * @param sequence String representing haplotype alignment.
 */
class Haplotype(val sequence: String, region: Seq[RichADAMRecord], val reference: String = "") {

  val hmm = new HMMAligner()
  val hasVariants = hmm.alignSequences(reference, sequence, null)
  val alignment = hmm.getAlignment()

  val perReadLikelihoods: Seq[Double] = region.map(read => {
    try {
      hmm.alignSequences(sequence, read.getSequence.toString, null)
      val readLikelihood = hmm.getLikelihood + hmm.getPrior
      readLikelihood
    } catch {
      case _: Throwable => {
        0.0
      }
    }
  })

  val readsLikelihood = perReadLikelihoods.sum

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
class HaplotypePair(val haplotype1: Haplotype, val haplotype2: Haplotype, hmm: HMMAligner) {

  val hasVariants = haplotype1.hasVariants || haplotype2.hasVariants
  val pairLikelihood = scorePairLikelihood

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
    hmm.alignSequences(haplotype2.sequence, haplotype1.sequence, null)
    val priorProb = hmm.getPrior
    readsProb + priorProb
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
    } else if (pair1.pairLikelihood > pair2.pairLikelihood) {
      1
    } else {
      0
    }
  }
}
