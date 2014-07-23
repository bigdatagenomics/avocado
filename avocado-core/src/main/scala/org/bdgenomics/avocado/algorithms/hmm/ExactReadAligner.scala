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

import org.bdgenomics.adam.algorithms.prefixtrie.DNAPrefixTrie
import org.bdgenomics.adam.util.PhredUtils
import scala.annotation.tailrec
import scala.math.log10

class ExactReadAligner(reference: String, rLen: Int) extends Aligner {

  // cut reference up
  protected val seq = reference.sliding(rLen).toArray
  protected val refLen = reference.length

  private def score(stat: (Double, Int), testSequence: String): Alignment = {
    val (alignScore, alignPos) = stat

    val matchSeq = seq(alignPos).zip(testSequence)
      .map(p => {
        if (p._1 == p._2) {
          "M"
        } else {
          "X"
        }
      }).reduce(_ + _)

    val alignedSequence = ("_" * alignPos) + testSequence + ("_" * (refLen - alignPos - rLen))
    val stateSequence = ("P" * alignPos) + matchSeq + ("P" * (refLen - alignPos - rLen))

    new Alignment(alignScore, 0.0, reference, alignedSequence, stateSequence)
  }

  /**
   * Aligns sequences.
   *
   * @param refSequence Reference sequence over the active region.
   * @param testSequence Sequence being scored.
   * @param testQualities String of qualities. Not currently used.
   * @return Alignment which stores the aligned sequences and likelihoods
   */
  def alignSequences(refSequence: String, testSequence: String, testQualities: String): Alignment = {
    val successQuals = testQualities.map(q => log10(PhredUtils.phredToSuccessProbability(q.toInt - 33)))
    val errorQuals = testQualities.map(q => log10(PhredUtils.phredToErrorProbability(q.toInt - 33)))

    val quals = seq.map(s => {
      (0 until rLen).map(i => {
        if (s(i) == testSequence(i)) {
          successQuals(i)
        } else {
          errorQuals(i)
        }
      }).reduce(_ + _)
    }).zipWithIndex
      .reduce((p1: (Double, Int), p2: (Double, Int)) => {
        if (p1._1 >= p2._1) {
          p1
        } else {
          p2
        }
      })

    score(quals, testSequence)
  }
}
