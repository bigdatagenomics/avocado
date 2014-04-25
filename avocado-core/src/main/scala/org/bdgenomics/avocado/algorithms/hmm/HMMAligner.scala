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

import scala.annotation.tailrec

object HMMAligner {
  val debug = false

  def align(refSequence: String, testSequence: String, testQualities: String): Alignment = {
    val hmm = new HMMAligner
    hmm.alignSequences(refSequence, testSequence, testQualities)
  }

}

/**
 * Pairwise alignment HMM. See the Durbin textbook (1998), chapter 4.
 *
 * All likelihoods are computed in log-space as are the input parameters
 */
class HMMAligner(val LOG_GAP_OPEN_PENALTY: Double = -4.0,
                 val LOG_GAP_CONTINUE_PENALTY: Double = -2.0,
                 val LOG_SNP_RATE: Double = -3.0,
                 val LOG_INDEL_RATE: Double = -4.0) {

  val transitionMatrix = new TransitionMatrix(LOG_GAP_OPEN_PENALTY,
    LOG_GAP_CONTINUE_PENALTY,
    LOG_SNP_RATE,
    LOG_INDEL_RATE)
  /**
   * Aligns sequences.
   *
   * @param refSequence Reference sequence over the active region.
   * @param testSequence Sequence being scored.
   * @param testQualities String of qualities. Not currently used.
   * @return Alignment which stores the aligned sequences and likelihoods
   */
  def alignSequences(refSequence: String, testSequence: String, testQualities: String): Alignment = {

    computePathLikelihood(refSequence, testSequence, testQualities)
    constructAlignment(refSequence, testSequence, transitionMatrix)
  }

  /**
   * Builds the aligned sequences by choosing the most likely state transitions
   *
   * @param refSequence Reference sequence over the active region.
   * @param testSequence Sequence being scored.
   * @return Alignment which stores the aligned sequences and likelihoods
   */
  def constructAlignment(refSequence: String, testSequence: String, transitionMatrix: TransitionMatrix): Alignment = {

    val alignmentLikelihood = transitionMatrix.getAlignmentLikelihood
    val paddedRefLen = refSequence.length + 1
    val paddedTestLen = testSequence.length + 1
    val stride = paddedRefLen

    @tailrec
    def constructAlignmentSequences(i: Int, j: Int, revAlignedRefSeq: String = "", revAlignedTestSeq: String = "", revAlignment: String = "", numSnps: Int = 0, numIndels: Int = 0): Alignment = {
      val idx = i * stride + j
      if (i <= 0 || j <= 0) {
        // Compute the prior probability of the alignments, with the Dindel numbers.
        val hasVariants: Boolean = numSnps > 0 || numIndels > 0
        val alignmentPrior = LOG_SNP_RATE * numSnps + LOG_INDEL_RATE * numIndels

        new Alignment(alignmentLikelihood, alignmentPrior, revAlignedRefSeq.reverse, revAlignedTestSeq.reverse, revAlignment.reverse, hasVariants)
      } else {
        /*
         * Check for scoring at each position, and then suggest next move. At this step, we identify
         * the next direction to move by looking at the "state" of the current coordinate. We call
         * our current 'state' by choosing the state with the highest cumulative likelihood.
         */
        val nextDirection = transitionMatrix.getMostLikelyState(idx)
        nextDirection match {
          case AlignmentState.Match => {
            val isSnp = testSequence(i - 1) != refSequence(j - 1)
            val alignmentCharacter = if (isSnp) 'X' else '='
            val addSnp = if (isSnp) 1 else 0
            constructAlignmentSequences(i - 1, j - 1, revAlignedRefSeq + refSequence(j - 1), revAlignedTestSeq + testSequence(i - 1), revAlignment + alignmentCharacter, numSnps + addSnp, numIndels)
          }
          case AlignmentState.Insertion => {
            // Inserted base from reference, add gap in reference alignment, check next base in test
            constructAlignmentSequences(i - 1, j, revAlignedRefSeq + '_', revAlignedTestSeq + testSequence(i - 1), revAlignment + 'I', numSnps, numIndels + 1)
          }
          case AlignmentState.Deletion => {
            // Deleted base from reference, add gap in text alignment, check next base in reference
            constructAlignmentSequences(i, j - 1, revAlignedRefSeq + refSequence(j - 1), revAlignedTestSeq + '_', revAlignment + 'D', numSnps, numIndels + 1)
          }
        }
      }

    }

    // Construct alignment sequence recursively starting at the end of the sequences
    constructAlignmentSequences(paddedTestLen - 1, paddedRefLen - 1)
  }

  /**
   * Compute the optimal transition matrix based on gap penalties
   *
   * @param refSequence Reference sequence over the active region.
   * @param testSequence Sequence being scored.
   * @param testQualities String of qualities. Not currently used.
   * @return TransitionMatrix which stores the likelihoods of each state and every pair of positions in the reference and test
   */
  private def computePathLikelihood(refSequence: String, testSequence: String, testQualities: String): TransitionMatrix = {
    val paddedRefLen = refSequence.length + 1
    val paddedTestLen = testSequence.length + 1
    val stride = paddedRefLen

    transitionMatrix.reallocate(paddedRefLen, paddedTestLen)

    for (testSeqPos <- 0 until paddedTestLen) {
      for (refSeqPos <- 0 until paddedRefLen) {
        val m: Double = transitionMatrix.getMatchLikelihood(testSeqPos, refSeqPos, stride, refSequence, testSequence)
        val ins: Double = transitionMatrix.getInsertionLikelihood(testSeqPos, refSeqPos, stride)
        val del: Double = transitionMatrix.getDeletionLikelihood(testSeqPos, refSeqPos, stride)
        if (testSeqPos > 0 || refSeqPos > 0) {
          val idx = testSeqPos * stride + refSeqPos
          transitionMatrix.matches(idx) = m
          transitionMatrix.inserts(idx) = ins
          transitionMatrix.deletes(idx) = del
        }
      }
    }
    transitionMatrix
  }
}
