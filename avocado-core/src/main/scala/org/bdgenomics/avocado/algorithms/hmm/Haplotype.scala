package org.bdgenomics.avocado.algorithms.hmm

import scala.collection.mutable.ArrayBuffer
import org.bdgenomics.adam.avro.ADAMRecord
import scala.math._
import scala.math.Ordering
/**
 * Haplotype generated from HMM alignment.
 *
 * @param sequence String representing haplotype alignment.
 */
class Haplotype(val sequence: String) {
  var perReadLikelihoods = new ArrayBuffer[Double]
  var readsLikelihood = Double.NegativeInfinity
  var hasVariants = false
  var alignment = new ArrayBuffer[(Int, Char)]

  /**
   * Score likelihood of reads when assembled into haplotype.
   *
   * @param hmm HMM aligner to use.
   * @param reads Sequence of reads to use.
   * @return Likelihood that reads are properly aligned.
   */
  def scoreReadsLikelihood(hmm: HMMAligner, reads: Seq[ADAMRecord]): Double = {
    perReadLikelihoods.clear
    readsLikelihood = 0.0
    for (r <- reads) {
      try {
        if (HMMAligner.debug) println(r.getSequence.toString + ", " + sequence)
        hmm.alignSequences(sequence, r.getSequence.toString, null)
        val readLike = hmm.getLikelihood + hmm.getPrior
        perReadLikelihoods += readLike
        readsLikelihood += readLike
      } catch {
        case _: Throwable => {
          perReadLikelihoods += 0.0
          readsLikelihood += 0.0
        }
      }
    }
    readsLikelihood
  }

  /**
   * Aligns reads to reference, and stores cigar details.
   *
   * @param hmm HMM aligner to use.
   * @param refHaplotype Haplotype for reference in this location.
   * @return True if region has variants.
   */
  def alignToReference(hmm: HMMAligner, refHaplotype: Haplotype): Boolean = {
    // TODO(peter, 12/8) store the alignment details (including the cigar).
    hasVariants = hmm.alignSequences(refHaplotype.sequence, sequence, null)
    alignment = hmm.getAlignment
    hasVariants
  }

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
   * likelihood than h2.
   *
   * @param h1 First haplotype to compare.
   * @param h2 Second haplotype to compare.
   * @return Ordering info for haplotypes.
   */
  def compare(h1: Haplotype, h2: Haplotype): Int = {
    if (h1.readsLikelihood < h2.readsLikelihood) {
      -1
    } else if (h1.readsLikelihood > h2.readsLikelihood) {
      1
    } else {
      0
    }
  }
}

/**
 * Class for a pairing of two haplotypes.
 *
 * @param haplotype1 First haplotype of pair.
 * @param haplotype2 Second haplotype of pair.
 */
class HaplotypePair(val haplotype1: Haplotype, val haplotype2: Haplotype) {
  var pairLikelihood = Double.NegativeInfinity
  var hasVariants = false

  override def toString(): String = {
    haplotype1.sequence + ", " + haplotype2.sequence + ", " + ("%1.3f" format pairLikelihood)
  }

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

  /**
   * Scores likelihood of two paired haplotypes and their alignment.
   *
   * @param hmm HMM aligner to use.
   * @param reads Sequence of reads that are evidence for haplotype.
   * @return Phred scaled likelihood.
   */
  def scorePairLikelihood(hmm: HMMAligner, reads: Seq[ADAMRecord]): Double = {
    var readsProb = 0.0
    for (i <- 0 until reads.length) {
      val readLikelihood1 = haplotype1.perReadLikelihoods(i)
      val readLikelihood2 = haplotype2.perReadLikelihoods(i)
      readsProb += exactLogSumExp10(readLikelihood1, readLikelihood2) - log10(2.0)
    }
    hmm.alignSequences(haplotype2.sequence, haplotype1.sequence, null)
    val priorProb = hmm.getPrior
    pairLikelihood = readsProb + priorProb
    pairLikelihood
  }

  /**
   * Aligns haplotype pair to reference. Joins variants of alignments.
   *
   * @param hmm HMM aligner to use.
   * @param refHaplotype Reference haplotype for active region.
   * @return True if region has variants.
   */
  def alignToReference(hmm: HMMAligner, refHaplotype: Haplotype): Boolean = {
    hasVariants = haplotype1.alignToReference(hmm, refHaplotype)
    if (haplotype2 != haplotype1) {
      hasVariants |= haplotype2.alignToReference(hmm, refHaplotype)
    }
    hasVariants
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
