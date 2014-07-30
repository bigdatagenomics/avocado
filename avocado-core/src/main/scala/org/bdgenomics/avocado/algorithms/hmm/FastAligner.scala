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

class FastAligner(reference: String, rLen: Int) extends Aligner {

  // cut into a trie
  protected val trie = DNAPrefixTrie(reference.sliding(rLen)
    .zipWithIndex
    .toMap)
  protected val refLen = reference.length
  protected val matchSeq = "M" * rLen

  private def score(stat: (Option[Int], Int), testQualities: String, testSequence: String): Alignment = {
    val (mismatchPos, alignPos) = stat

    val quals = testQualities.toArray
      .zipWithIndex

    val (eQual, mQualArray, stateSeq) = if (mismatchPos.isDefined) {
      val mp = mismatchPos.get
      (log10(PhredUtils.phredToErrorProbability(quals(mp)._1.toInt - 33)),
        quals.filter(kv => kv._2 != mp),
        ("M" * mp) + "X" + ("M" * (rLen - mp - 1)))
    } else {
      (0.0, quals, matchSeq)
    }

    val qual = eQual + mQualArray.map(c => log10(PhredUtils.phredToSuccessProbability(c._1.toInt - 33)))
      .reduce(_ + _)

    val alignedSequence = ("_" * alignPos) + testSequence + ("_" * (refLen - alignPos - rLen))
    val stateSequence = ("P" * alignPos) + stateSeq + ("P" * (refLen - alignPos - rLen))

    new Alignment(qual, 0.0, reference, alignedSequence, stateSequence)
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
    @tailrec def missedAlignmentHelper(pos: Iterator[Int], testSeq: String): (Option[Int], Int) = {
      if (!pos.hasNext) {
        throw new IllegalArgumentException("Couldn't find alignment")
      } else {
        val p = pos.next
        val entries = trie.search(testSeq.take(p) + "*" + testSeq.drop(p + 1))

        if (!entries.isEmpty) {
          (Some(p), entries.head._2)
        } else {
          missedAlignmentHelper(pos, testSeq)
        }
      }
    }

    val alignmentStat = trie.getIfExists(testSequence)
      .fold(missedAlignmentHelper(testQualities.toSeq
        .zipWithIndex
        .sortBy(kv => kv._1)
        .map(kv => kv._2)
        .toIterator, testSequence))(p => (None.asInstanceOf[Option[Int]], p))

    score(alignmentStat, testQualities, testSequence)
  }
}
