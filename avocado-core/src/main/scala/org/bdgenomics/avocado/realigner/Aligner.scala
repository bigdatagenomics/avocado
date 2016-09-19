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
package org.bdgenomics.avocado.realigner

import org.bdgenomics.avocado.models.{
  Deletion,
  Insertion,
  Match,
  ObservationOperator
}
import scala.annotation.tailrec

/**
 * Core code for aliging two sequences (read and ref) against each other.
 */
object Aligner {

  /**
   * Trims matching bases off of two sequences that contain an INDEL.
   *
   * @param ref First sequence for comparison.
   * @param alt Second sequence for comparison.
   * @return Returns a tuple containing (trimmed ref, trimmed alt, number of
   *   bases trimmed from start, number of bases trimmed from end).
   */
  private[realigner] def zipAndTrim(ref: String,
                                    alt: String): (String, String, Int, Int) = {

    def doZipAndTrim(refSeq: String, altSeq: String): (String, String, Int) = {
      val zippedSequences = refSeq.zip(altSeq)

      // we count off the number of bases at the start/end that match for trimming
      val matchIdx = zippedSequences.indexWhere(p => p._1 != p._2)

      // if we've already trimmed once, we may not find an index on the second
      // trimming pass and matchIdx will be -1
      val startMatch = if (matchIdx == -1) {
        zippedSequences.length
      } else {
        matchIdx
      }

      // trim and return
      (refSeq.drop(startMatch),
        altSeq.drop(startMatch),
        startMatch)
    }

    // reverse and zip and trim
    val (revRef, revAlt, endTrimmed) = doZipAndTrim(ref.reverse, alt.reverse)

    // reverse and zip and trim forward
    val (fwdRef, fwdAlt, startTrimmed) = doZipAndTrim(revRef.reverse, revAlt.reverse)

    (fwdRef, fwdAlt, startTrimmed, endTrimmed)
  }

  /**
   * @param sequence Sequence to cut into k-mers.
   * @param kmerLength The length k of the k-mers.
   * @return Returns a map between
   */
  private[realigner] def toKmers(sequence: String,
                                 kmerLength: Int): Map[String, Int] = {

    // we can't cut up a sequence that is shorter than the k-mer length,
    // so don't even try!
    if (sequence.length < kmerLength) {
      Map.empty
    } else {
      // cut up the k-mers and attach their location in the sequence
      val kmerMap = sequence.sliding(kmerLength)
        .zipWithIndex
        .toMap

      assert((sequence.size - kmerLength + 1) == kmerMap.size,
        "Input sequence contains a repeat.")
      kmerMap
    }
  }

  /**
   * Aligns two sequences that have been trimmed and partitioned.
   *
   * @param trimmedRef The trimmed/partitioned reference sequence.
   * @param trimmedAlt The trimmed/partitioned alt sequence.
   * @return Returns the local alignment of these two sequence segments.
   */
  private def alignSegment(trimmedRef: String,
                           trimmedAlt: String): Iterable[ObservationOperator] = {

    val refLength = trimmedRef.length
    val altLength = trimmedAlt.length

    // simple snp or mnp
    if (refLength == altLength) {

      // per base in the snp/mnp, check equality
      val verboseOperators = trimmedRef.zip(trimmedAlt)
        .map(p => {
          val (r, a) = p

          if (r == a) {
            Match(1)
          } else {
            Match(1, Some(r.toString))
          }
        })

      // collapse down the operators
      ObservationOperator.collapse(verboseOperators)
    } else {
      // insertion or deletion

      // three possible cases:
      //
      // - ref is empty --> simple insertion
      // - alt is empty --> simple deletion
      // - both have a value --> complex indel
      if (refLength == 0) {
        Iterable(Insertion(altLength))
      } else if (altLength == 0) {
        Iterable(Deletion(trimmedRef))
      } else {
        // complex indel = both insertion and deletion (or indel + mnp)
        // put longer event "first"
        if (altLength > refLength) {
          Iterable(Insertion(altLength),
            Deletion(trimmedRef))
        } else {
          Iterable(Deletion(trimmedRef),
            Insertion(altLength))
        }
      }
    }
  }

  /**
   * Aligns two sequences against each other.
   *
   * Performs trimming and partitioning before aligning the two sequences.
   *
   * @param ref Reference sequence to align against.
   * @param alt Alternate sequence to align against ref.
   * @param kmerLength Length k to use for cutting up k-mers.
   * @return Returns a sequence of alignment operators describing the pairwise
   *   alignment of the two sequences.
   */
  def align(ref: String,
            alt: String,
            kmerLength: Int): Seq[ObservationOperator] = {

    // trim ref and alt
    // drop the first and last matching bases from the sequences
    val (trimmedRef, trimmedAlt, startLength, endLength) = zipAndTrim(ref, alt)

    // this yields the matches at the start and the end
    val startEvent = Match(startLength)
    val endEvent = Match(endLength)

    // make k-mer maps
    val refKmerMap = toKmers(trimmedRef, kmerLength)
    val altKmerMap = toKmers(trimmedAlt, kmerLength)

    // find k-mer set intersection
    val intersect = refKmerMap.keySet & altKmerMap.keySet

    def emit: Seq[ObservationOperator] = {
      // align the segment we have and emit start/end events
      Seq.concat(Some(startEvent),
        alignSegment(trimmedRef, trimmedAlt),
        Some(endEvent))
    }

    // if we find no intersecting k-mers, then we are done!
    if (intersect.isEmpty) {

      emit
    } else {

      // extract indices
      val indices = kmerIndices(intersect,
        refKmerMap,
        altKmerMap)

      // and validate the concordance of the indices we got back
      // if we have a non-concordant set of indices, we can't clean up any
      // further and will just emit our alignment
      if (indicesHaveConcordantOrder(indices)) {

        // get the blocks that our ranges represent
        val ranges = indicesToBlocks(indices,
          trimmedRef,
          trimmedAlt,
          kmerLength)

        // convert these blocks into alignments by folding
        val alignmentOperators = ranges.flatMap(_.fold(alignSegment))

        Seq.concat(Some(startEvent),
          alignmentOperators,
          Some(endEvent))
      } else {
        emit
      }
    }
  }

  /**
   * @param indices Seq of (ref index, alt index) tuples for k-mers that are in
   *   both ref and alt sequences.
   * @param refSeq The reference sequence.
   * @param altSeq The alternate sequence.
   * @param kmerLength The length k of the k-mers.
   * @return Returns a seq of blocks that can be folded over to complete the
   *   alignment.
   *
   * @note indices is assumed to be sorted and concordant.
   */
  private[realigner] def indicesToBlocks(indices: Seq[(Int, Int)],
                                         refSeq: String,
                                         altSeq: String,
                                         kmerLength: Int): Seq[Block] = {

    @tailrec def chopBlocks(indexIter: Iterator[(Int, Int)],
                            ref: String,
                            alt: String,
                            matchLength: Int = -1,
                            blockList: List[Block] = List.empty): Seq[Block] = {

      if (!indexIter.hasNext) {

        val blocks = if (matchLength > 0) {
          MatchBlock(matchLength) :: blockList
        } else {
          blockList
        }

        (UnknownBlock(ref, alt) :: blocks).toSeq
      } else {

        // pop from the iterator
        val (nextRefIdx, nextAltIdx) = indexIter.next

        // was our last match > 0? if so, prepend
        def blocks: List[Block] = if (matchLength > 0) {
          MatchBlock(matchLength) :: blockList
        } else {
          blockList
        }

        // our sequences should either be:
        // 1. len > idx + kmerLength, or
        // 2. len == idx
        val refExtension = ref.length - nextRefIdx
        val altExtension = alt.length - nextAltIdx

        val nextRef = ref.take(nextRefIdx)
        val nextAlt = alt.take(nextAltIdx)
        val (newMatchLength, newList) = if (refExtension >= kmerLength && altExtension >= kmerLength) {
          (kmerLength,
            UnknownBlock(ref.takeRight(refExtension - kmerLength),
              alt.takeRight(altExtension - kmerLength)) :: blocks)
        } else if (refExtension > 1 && altExtension >= 1 &&
          refExtension > altExtension) {
          (altExtension,
            (UnknownBlock(ref.takeRight(refExtension - altExtension), "")
              :: blocks))
        } else if (altExtension > 1 && refExtension >= 1 &&
          altExtension > refExtension) {
          (refExtension,
            (UnknownBlock("", alt.takeRight(altExtension - refExtension))
              :: blocks))
        } else {
          assert(refExtension == 1 && altExtension == 1)
          (matchLength + 1,
            blockList)
        }

        chopBlocks(indexIter,
          nextRef, nextAlt,
          matchLength = newMatchLength,
          blockList = newList)
      }
    }

    chopBlocks(indices.reverseIterator,
      refSeq,
      altSeq)
  }

  /**
   * @param intersection The set of k-mers seen in both the ref and the alt.
   * @param refMap A map from ref k-mers to their sequence index.
   * @param altMap A map from alt k-mers to their sequence index.
   * @return Returns a sorted seq containing tuples of the sequence indices for
   *   each k-mer in both the ref and the alt.
   */
  private[realigner] def kmerIndices(intersection: Set[String],
                                     refMap: Map[String, Int],
                                     altMap: Map[String, Int]): Seq[(Int, Int)] = {

    intersection.map(k => {
      (refMap(k), altMap(k))
    }).toSeq
      .sortBy(_._1)
  }

  /**
   * Two index sets are concordant if the indices have the same order in each.
   *
   * @param indices A sorted seq containing tuples of indices.
   * @return Returns true if after sorting the set on the first value in the
   *   index tuple, the second value is ordered as well.
   */
  private[realigner] def indicesHaveConcordantOrder(
    indices: Seq[(Int, Int)]): Boolean = {

    @tailrec def isSecondarilySorted(iter: Iterator[(Int, Int)],
                                     lastIdx: Int = -1): Boolean = {
      if (!iter.hasNext) {

        // if we have made it to the end, then our iterator was monotonic
        true
      } else {

        // pop the index from the head of the list
        val (_, idx) = iter.next

        // by definition, indices must be >= 0
        assert(idx >= 0)

        // if our index is smaller than our last index, then the list is not
        // increasing monotonically. if it is, then keep recursing.
        if (idx <= lastIdx) {
          false
        } else {
          isSecondarilySorted(iter, idx)
        }
      }
    }

    // check for second index sort
    isSecondarilySorted(indices.toIterator)
  }
}
