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
import org.bdgenomics.avocado.algorithms.hmm.AlignmentState._

class TransitionMatrix(val LOG_GAP_OPEN: Double = -4.0,
                       val LOG_GAP_CONTINUE: Double = -2.0,
                       val LOG_SNP_RATE: Double = -3.0,
                       val LOG_INDEL_RATE: Double = -4.0) {

  var matSize = 0
  var stride = 0
  var matches = new Array[Double](matSize)
  var inserts = new Array[Double](matSize)
  var deletes = new Array[Double](matSize)

  // P( match -> match) = 1 - 2 * P(match -> indel) = 1 - 2 * P( gap-open )
  // in log-space
  val LOG_MATCH_TO_MATCH = log10(1 - 2 * pow(10, LOG_GAP_OPEN))

  // P( indel -> match) = 1 - P(indel -> indel) = 1 - P( gap-continue)
  // in log-space
  val LOG_GAP_CLOSE = log10(1 - pow(10, LOG_GAP_CONTINUE))

  // P( match) = 1 - P(mismatch)
  val MATCH_PROB = 1 - pow(10, LOG_SNP_RATE)

  // in log-space
  val LOG_MATCH_PROB = log10(MATCH_PROB)

  /**
   *
   *  This creates a single backing transition matrix and only reallocates
   *  the matrices if the length of the sequences has increased
   *
   */
  def reallocate(rows: Int, columns: Int) = {
    val newMatrixSize = rows * columns
    stride = rows
    if (newMatrixSize > matSize) {
      matches = new Array[Double](newMatrixSize)
      inserts = new Array[Double](newMatrixSize)
      deletes = new Array[Double](newMatrixSize)
    }

    matSize = newMatrixSize
    initialise(rows)
  }

  def getAlignmentLikelihood(): Double = {
    max(matches(matSize - 1), max(inserts(matSize - 1), deletes(matSize - 1)))
  }

  def initialise(sequenceLength: Int) = {
    matches(0) = 2 * log10(1.0 + sequenceLength)
    inserts(0) = Double.NegativeInfinity
    deletes(0) = Double.NegativeInfinity
  }

  def getMostLikelyState(idx: Int, epsilon: Double = 1e-2): AlignmentState = {
    if (matches(idx) + epsilon >= inserts(idx) && matches(idx) + epsilon >= deletes(idx))
      AlignmentState.Match
    else if (inserts(idx) >= deletes(idx)) {
      AlignmentState.Insertion
    } else {
      AlignmentState.Deletion
    }
  }

  def getInsertionLikelihood(i: Int, j: Int, stride: Int): Double = {
    if (i >= 1) {
      val idx = (i - 1) * stride + j
      max(gapProbability(AlignmentState.Insertion, idx), gapProbability(AlignmentState.Match, idx))
    } else {
      Double.NegativeInfinity
    }
  }

  def getDeletionLikelihood(i: Int, j: Int, stride: Int): Double = {
    if (j >= 1) {
      val idx = i * stride + (j - 1)
      max(gapProbability(AlignmentState.Deletion, idx), gapProbability(AlignmentState.Match, idx))
    } else {
      Double.NegativeInfinity
    }
  }

  def getMatchLikelihood(i: Int, j: Int, stride: Int, refSequence: String, testSequence: String): Double = {
    if (i >= 1 && j >= 1) {
      val testBase = testSequence(i - 1)
      val refBase = refSequence(j - 1)
      val (emitProbability) = if (testBase == refBase) {
        LOG_MATCH_PROB
      } else {
        LOG_SNP_RATE
      }
      val idx = (i - 1) * stride + (j - 1)
      val p = AlignmentState.values.map(state => matchProbability(state, idx)).max + emitProbability
      p
    } else {
      Double.NegativeInfinity
    }
  }

  private def gapProbability(state: AlignmentState, idx: Int): Double = {
    state match {
      case AlignmentState.Insertion                       => inserts(idx) + LOG_GAP_CONTINUE
      case AlignmentState.Deletion                        => deletes(idx) + LOG_GAP_CONTINUE
      case AlignmentState.Match | AlignmentState.Mismatch => matches(idx) + LOG_GAP_OPEN
    }
  }

  private def matchProbability(state: AlignmentState, idx: Int): Double = {
    state match {
      case AlignmentState.Insertion                       => inserts(idx) + LOG_GAP_CLOSE
      case AlignmentState.Deletion                        => deletes(idx) + LOG_GAP_CLOSE
      case AlignmentState.Match | AlignmentState.Mismatch => matches(idx) + LOG_MATCH_TO_MATCH
    }
  }
}
