/*
 * Copyright (c) 2013-2014. Regents of the University of California
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

import org.bdgenomics.adam.avro.ADAMRecord
import scala.collection.mutable.ArrayBuffer
import scala.math._

object HMMAligner {
  val debug = false
}

/**
 * Pairwise alignment HMM. See the Durbin textbook (1998), chapter 4.
 */
class HMMAligner {
  var refSequence: String = null
  var refAligned: String = null

  var testSequence: String = null
  var testAligned: String = null

  var testQualities: String = null

  var refOffset = 0

  var paddedRefLen = 0
  var paddedTestLen = 0

  var stride = 0
  var matSize = 1

  var eta = Double.NegativeInfinity

  // This uses the quick and dirty numbers from the Dindel (2011) paper.
  // TODO(peter, 12/7) I'm forgetting a factor of 3 somewhere...
  val mismatchPrior = -3.0 - log10(3.0)
  val matchPrior = log10(1.0 - 1.0e-3)
  val indelPrior = -4.0
  val indelToMatchPrior = -4.0
  val indelToIndelPrior = -4.0

  private var matches: Array[Double] = Array(0.0)
  private var inserts: Array[Double] = Array(0.0)
  private var deletes: Array[Double] = Array(0.0)

  private var traceMatches: Array[Char] = Array('\0')
  private var traceInserts: Array[Char] = Array('\0')
  private var traceDeletes: Array[Char] = Array('\0')

  private var alignment: ArrayBuffer[(Int, Char)] = null

  private var alignmentLikelihood = Double.NegativeInfinity
  private var alignmentPrior = Double.NegativeInfinity

  /**
   * Aligns sequences.
   *
   * @param refSequence Reference sequence over the active region.
   * @param testSequence Sequence being scored.
   * @param testQualities String of qualities. Not currently used.
   * @return True if sequence has variants.
   */
  def alignSequences(refSequence: String, testSequence: String, testQualities: String): Boolean = {

    def fpCompare(a: Double, b: Double, eps: Double): Boolean = abs(a - b) <= eps

    paddedRefLen = refSequence.length + 1
    paddedTestLen = testSequence.length + 1
    stride = paddedRefLen
    val oldMatSize = matSize
    matSize = max(matSize, paddedTestLen * paddedRefLen)
    if (matSize > oldMatSize && matSize > 0) {
      matches = new Array[Double](matSize)
      inserts = new Array[Double](matSize)
      deletes = new Array[Double](matSize)
      traceMatches = new Array[Char](matSize)
      traceInserts = new Array[Char](matSize)
      traceDeletes = new Array[Char](matSize)
    }
    else if (matSize <= 1 || oldMatSize <= 1) {
      matches = new Array[Double](1)
      inserts = new Array[Double](1)
      deletes = new Array[Double](1)
      traceMatches = new Array[Char](1)
      traceInserts = new Array[Char](1)
      traceDeletes = new Array[Char](1)
    }

    // Note: want to use the _test_ haplotype length here, not the ref length.
    eta = -log10(1.0 + testSequence.length)

    // Compute the optimal alignment.
    // TODO(peter, 12/4) shortcut b/w global and local alignment: use a custom
    // start position in the reference haplotype.
    matches(0) = 2.0 * eta
    inserts(0) = Double.NegativeInfinity
    deletes(0) = Double.NegativeInfinity
    for (i <- 0 until paddedTestLen) {
      for (j <- 0 until paddedRefLen) {
        if (i > 0 || j > 0) {
          val (m, trM) = if (i >= 1 && j >= 1) {
            val testBase = testSequence(i - 1)
            val refBase = refSequence(j - 1)
            // TODO(peter, 12/7) there is a second constant term to the prior...
            val prior = if (testBase == refBase) {
              matchPrior / (indelPrior * indelPrior)
            }
            else {
              mismatchPrior / (indelPrior * indelPrior)
            }
            val idx = (i - 1) * stride + (j - 1)
            val mMatch = matches(idx)
            val mInsert = inserts(idx)
            val mDelete = deletes(idx)
            val mMax = max(mMatch, max(mInsert, mDelete)) + prior
            val t = if (fpCompare(mMax, mMatch, 0.01)) {
              'M'
            }
            else if (fpCompare(mMax, mInsert, 0.01)) {
              'I'
            }
            else if (fpCompare(mMax, mDelete, 0.01)) {
              'D'
            }
            else {
              '.'
            }
            (mMax, t)
          }
          else {
            (Double.NegativeInfinity, '.')
          }
          val (ins, trIns) = if (i >= 1) {
            val idx = (i - 1) * stride + j
            val insMatch = matches(idx) + indelToMatchPrior
            val insInsert = inserts(idx) + indelToIndelPrior
            val insMax = max(insMatch, insInsert)
            val t = if (fpCompare(insMax, insMatch, 0.01)) {
              'M'
            }
            else if (fpCompare(insMax, insInsert, 0.01)) {
              'I'
            }
            else {
              '.'
            }
            (insMax, t)
          }
          else {
            (Double.NegativeInfinity, '.')
          }
          val (del, trDel) = if (j >= 1) {
            val idx = i * stride + (j - 1)
            val delMatch = matches(idx) + indelToMatchPrior
            val delDelete = deletes(idx) + indelToIndelPrior
            val delMax = max(delMatch, delDelete)
            val t = if (fpCompare(delMax, delMatch, 0.01)) {
              'M'
            }
            else if (fpCompare(delMax, delDelete, 0.01)) {
              'D'
            }
            else {
              '.'
            }
            (delMax, t)
          }
          else {
            (Double.NegativeInfinity, '.')
          }
          val idx = i * stride + j
          matches(idx) = m
          inserts(idx) = ins
          deletes(idx) = del
          traceMatches(idx) = trM
          traceInserts(idx) = trIns
          traceDeletes(idx) = trDel
        }
      }
    }

    if (matSize > 0) {
      alignmentLikelihood = max(matches(matSize - 1), max(inserts(matSize - 1), deletes(matSize - 1)))
    }
    else {
      alignmentLikelihood = max(matches(0), max(inserts(0), deletes(0)))
    }

    def printArray(a: Array[Double]) {
      for (i <- 0 until paddedTestLen) {
        var s = ""
        for (j <- 0 until paddedRefLen) {
          val idx = i * stride + j
          s += "%2.2f" format a(idx)
          s += "\t"
        }
        if (HMMAligner.debug) println(s)
      }
    }

    if (HMMAligner.debug) {
      println("matches")
      printArray(matches)
      println("inserts")
      printArray(inserts)
      println("deletes")
      printArray(deletes)
    }

    // Traceback to get the aligned sequences.
    var hasVariants = false
    var numSnps = 0
    var numIndels = 0
    var revAlignment = ""
    var revAlignedTestSeq = ""
    var revAlignedRefSeq = ""

    def indexMax(a: (Double, Int), b: (Double, Int)): (Double, Int) = {
      if (a._1 > b._1) {
        a
      }
      else {
        b
      }
    }

    def getArrayMax(array: Array[Double]): (Double, Int) = {
      array.zipWithIndex.reduce(indexMax)
    }

    val highestM = getArrayMax(matches)
    val highestI = getArrayMax(inserts)
    val highestD = getArrayMax(deletes)
    val highestIdx = indexMax(highestM, indexMax(highestI, highestD))._2

    var i = paddedTestLen - 1
    var j = paddedRefLen - 1

    if (HMMAligner.debug)
      println(i + ", " + j)

    while (i > 0 && j > 0) {
      val idx = i * stride + j

      if (HMMAligner.debug) {
        println(i + ", " + j + (": %1.1f" format matches(idx)) + (", %1.1f" format inserts(idx)) +
          (", %1.1f" format deletes(idx)))
      }

      val bestScore = max(matches(idx), max(inserts(idx), deletes(idx)))
      // TODO(peter, 12/7) here, build the aligned sequences.
      val tr = if (fpCompare(bestScore, matches(idx), 0.0001)) {
        revAlignedTestSeq += testSequence(i - 1)
        revAlignedRefSeq += refSequence(j - 1)
        if (testSequence(i - 1) != refSequence(j - 1)) {
          hasVariants = true
          numSnps += 1
          revAlignment += 'X'
        }
        else {
          revAlignment += '='
        }
        // FIXME
        traceMatches(idx)
      }
      else if (fpCompare(bestScore, inserts(idx), 0.0001)) {
        revAlignedTestSeq += testSequence(i - 1)
        revAlignedRefSeq += '_'
        hasVariants = true
        numIndels += 1
        revAlignment += 'I'
        traceInserts(idx)
      }
      else if (fpCompare(bestScore, deletes(idx), 0.0001)) {
        revAlignedRefSeq += refSequence(j - 1)
        revAlignedTestSeq += '_'
        hasVariants = true
        numIndels += 1
        revAlignment += 'D'
        traceDeletes(idx)
      }
      else {
        revAlignedTestSeq += 'x'
        revAlignedRefSeq += 'x'
        '.'
      }

      tr match {
        case 'M' => {
          i -= 1
          j -= 1
        }
        case 'I' => {
          i -= 1
        }
        case 'D' => {
          j -= 1
        }
        case _ => {
          i -= 1
          j -= 1
        } // TODO(peter, 12/8) the alignment is bad (probably a bug).
      }
    }

    testAligned = revAlignedTestSeq.reverse
    refAligned = revAlignedRefSeq.reverse

    if (HMMAligner.debug) {
      println("ta: " + testAligned)
      println("ra: " + refAligned)
    }

    var unitAlignment = revAlignment.reverse
    if (HMMAligner.debug) println(unitAlignment)
    alignment = new ArrayBuffer[(Int, Char)]
    var alignSpan: Int = 0
    var alignMove: Char = '.'
    for (i <- 0 until unitAlignment.length) {
      val move = unitAlignment(i)
      if (move != alignMove) {
        if (alignSpan > 0) {
          val tok = (alignSpan, alignMove)
          alignment += tok
        }
        alignSpan = 1
        alignMove = move
      }
      else {
        alignSpan += 1
      }
    }
    val tok = (alignSpan, alignMove)
    alignment += tok

    // Compute the prior probability of the alignments, with the Dindel numbers.
    alignmentPrior = mismatchPrior * numSnps + indelPrior * numIndels
    if (HMMAligner.debug) println("ap: " + alignmentPrior)

    // Returns whether the alignment has any SNPs or indels.
    hasVariants
  }

  /**
   * Compute the (log10) likelihood of aligning the test sequence to the ref.
   *
   * @return Log10 likelihood of alignment.
   */
  def getLikelihood(): Double = alignmentLikelihood

  /**
   * Compute the (log10) prior prob of observing the given alignment.
   *
   * @return Prior for alignment.
   */
  def getPrior(): Double = alignmentPrior

  /**
   * Compute the alignment tokens (equivalent to cigar).
   *
   * @return Alignment tokens.
   */
  def getAlignment(): ArrayBuffer[(Int, Char)] = {
    alignment.clone
  }

  /**
   * Return the aligned sequences.
   *
   * @return Tuple of sequences aligned to (ref, test) sequences.
   */
  def getAlignedSequences(): (String, String) = (refAligned, testAligned)
}

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
        val readLike = hmm.getLikelihood // - hmm.getPriora
        perReadLikelihoods += readLike
        readsLikelihood += readLike
      }
      catch {
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
    }
    else if (h1.readsLikelihood > h2.readsLikelihood) {
      1
    }
    else {
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
    }
    else if (pair1.pairLikelihood > pair2.pairLikelihood) {
      1
    }
    else {
      0
    }
  }
}
