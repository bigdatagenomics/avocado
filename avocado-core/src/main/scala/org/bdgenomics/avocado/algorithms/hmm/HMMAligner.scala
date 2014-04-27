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
object AlignmentState extends Enumeration {
  type AlignmentState = Value
  val Insertion, Match, Mismatch, Deletion, None = Value
}

/**
 * Pairwise alignment HMM. See the Durbin textbook (1998), chapter 4.
 */
class HMMAligner(val indelPrior: Double = -4.0,
                 val mismatchPrior: Double = -3.0 - log10(3.0),
                 val matchPrior: Double = log10(1.0 - 1.0e-3),
                 val indelToMatchPrior: Double = -4.0,
                 val indelToIndelPrior: Double = -4.0) {

  var refAligned: String = null

  var testAligned: String = null

  var matSize = 1

  private var matches: Array[Double] = Array(0.0)
  private var inserts: Array[Double] = Array(0.0)
  private var deletes: Array[Double] = Array(0.0)

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
    val paddedRefLen = refSequence.length + 1
    val paddedTestLen = testSequence.length + 1
    val stride = paddedRefLen
    computeAlignmentLikelihood(testSequence, refSequence)

    // Traceback to get the aligned sequences.
    var hasVariants = false
    var numSnps = 0
    var numIndels = 0
    var revAlignment = ""
    var revAlignedTestSeq = ""
    var revAlignedRefSeq = ""

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

      /*
       * Check for scoring at each position, and then suggest next move. At this step, we identify
       * the next direction to move by looking at the "state" of the current coordinate. We call
       * our current 'state' by choosing the state with the highest cumulative likelihood.
       */
      if (matches(idx) >= inserts(idx) && matches(idx) >= deletes(idx)) {
        revAlignedTestSeq += testSequence(i - 1)
        revAlignedRefSeq += refSequence(j - 1)
        if (testSequence(i - 1) != refSequence(j - 1)) {
          hasVariants = true
          numSnps += 1
          revAlignment += 'X'
        } else {
          revAlignment += '='
        }
        i -= 1
        j -= 1
      } else if (inserts(idx) >= deletes(idx)) {
        revAlignedTestSeq += testSequence(i - 1)
        revAlignedRefSeq += '_'
        hasVariants = true
        numIndels += 1
        revAlignment += 'I'
        i -= 1
      } else {
        revAlignedRefSeq += refSequence(j - 1)
        revAlignedTestSeq += '_'
        hasVariants = true
        numIndels += 1
        revAlignment += 'D'
        j -= 1
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
      } else {
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

  def computeAlignmentLikelihood(testSequence: String, refSequence: String) {
    val paddedRefLen = refSequence.length + 1
    val paddedTestLen = testSequence.length + 1
    val stride = paddedRefLen
    val oldMatSize = matSize
    matSize = max(matSize, paddedTestLen * paddedRefLen)

    if (matSize > oldMatSize) {
      matches = new Array[Double](matSize)
      inserts = new Array[Double](matSize)
      deletes = new Array[Double](matSize)
    }

    // Note: want to use the _test_ haplotype length here, not the ref length.
    val eta = -log10(1.0 + testSequence.length)

    // Compute the optimal alignment.
    // TODO(peter, 12/4) shortcut b/w global and local alignment: use a custom
    // start position in the reference haplotype.
    matches(0) = 2.0 * eta
    inserts(0) = Double.NegativeInfinity
    deletes(0) = Double.NegativeInfinity
    for (i <- 0 until paddedTestLen) {
      for (j <- 0 until paddedRefLen) {
        if (i > 0 || j > 0) {
          val m = if (i >= 1 && j >= 1) {
            val testBase = testSequence(i - 1)
            val refBase = refSequence(j - 1)
            // TODO(peter, 12/7) there is a second constant term to the prior...
            val prior = if (testBase == refBase) {
              matchPrior / (indelPrior * indelPrior)
            } else {
              mismatchPrior / (indelPrior * indelPrior)
            }
            val idx = (i - 1) * stride + (j - 1)
            val mMatch = matches(idx)
            val mInsert = inserts(idx)
            val mDelete = deletes(idx)
            val mMax = max(mMatch, max(mInsert, mDelete)) + prior
            mMax
          } else {
            Double.NegativeInfinity
          }
          val ins = if (i >= 1) {
            val idx = (i - 1) * stride + j
            val insMatch = matches(idx) + indelToMatchPrior
            val insInsert = inserts(idx) + indelToIndelPrior
            val insMax = max(insMatch, insInsert)
            insMax
          } else {
            Double.NegativeInfinity
          }
          val del = if (j >= 1) {
            val idx = i * stride + (j - 1)
            val delMatch = matches(idx) + indelToMatchPrior
            val delDelete = deletes(idx) + indelToIndelPrior
            val delMax = max(delMatch, delDelete)
            delMax
          } else {
            Double.NegativeInfinity
          }
          val idx = i * stride + j
          matches(idx) = m
          inserts(idx) = ins
          deletes(idx) = del
        }
      }
    }

    alignmentLikelihood = max(matches(matSize - 1), max(inserts(matSize - 1), deletes(matSize - 1)))
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
