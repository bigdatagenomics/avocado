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

import htsjdk.samtools.{ CigarOperator, TextCigarCodec }
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.util.MdTag
import org.bdgenomics.avocado.models.{
  Clipped,
  Deletion,
  Insertion,
  Match,
  ObservationOperator
}
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.utils.misc.Logging
import scala.annotation.tailrec
import scala.collection.JavaConversions._

/**
 * Singleton object for realigning reads against the reference genome.
 */
object Realigner extends Logging {

  /**
   * Realigns a set of reads against the reference genome.
   *
   * @param reads Reads to realign.
   * @param kmerLength The length k of the k-mers.
   * @return Returns the realigned reads.
   */
  def realign(reads: AlignmentRecordRDD,
              kmerLength: Int): AlignmentRecordRDD = {
    reads.transform(realignRdd(_, kmerLength))
  }

  /**
   * The realign method works on a wrapped genomic RDD. This is the actual
   * method that is used to transform the reads.
   *
   * This method maps over each read. On each read, we check to see if the read
   * is a candidate for realignment. If the read is a realignment candidate, we
   * go ahead and realign it. If not, we directly emit the input read.
   *
   * @param reads Reads to realign.
   * @param kmerLength The length k of the k-mers.
   * @return Returns the realigned reads.
   */
  private def realignRdd(reads: RDD[AlignmentRecord],
                         kmerLength: Int): RDD[AlignmentRecord] = {
    reads.map(r => {

      // we can't realign unmapped reads
      if (!r.getReadMapped) {
        r
      } else {

        // convert the cigar and MD tag into alignment operators
        val alignment = extractAlignmentOperators(r)

        // check whether the read should be realigned. if so, realign.
        if (alignment.isEmpty) {

          // if the alignment operator is empty, then we couldn't extract from the
          // cigar/md tag. log a warning message!
          log.warn(
            "Couldn't extract alignment from CIGAR (%s) and MD (%s) for %s. Skipping realignment.",
            r.getCigar,
            r.getMismatchingPositions,
            r.getReadName)

          r
        } else {

          // get the realignment blocks for this read
          val realignmentBlocks = RealignmentBlock(r.getSequence,
            alignment,
            kmerLength)

          // is this read good for realignment?
          if (isRealignmentCandidate(realignmentBlocks)) {
            realignRead(r, realignmentBlocks, kmerLength)
          } else {
            r
          }
        }
      }
    })
  }

  /**
   * Realigns a single read versus the reference.
   *
   * Extracts the reference, segments out the alignment blocks, and realigns.
   *
   * @param read Read to realign.
   * @param alignment The candidate realignment blocks.
   * @param kmerLength The length k of the k-mers.
   * @return Returns a read
   */
  private[realigner] def realignRead(
    read: AlignmentRecord,
    realignmentBlocks: Iterable[RealignmentBlock],
    kmerLength: Int): AlignmentRecord = {

    // realign blocks
    val realignedOperators = realignmentBlocks.flatMap(block => {

      def extractAndRealign(read: String,
                            ops: Iterable[ObservationOperator]): Iterable[ObservationOperator] = {

        // extract the reference sequence
        val ref = ObservationOperator.extractReference(read, ops)

        // align the read against the ref
        Aligner.align(ref, read, kmerLength)
      }

      block.fold(extractAndRealign)
    })

    // we may have a run of likewise blocks, so collapse on down
    val collapsedOperators = ObservationOperator.collapse(realignedOperators)

    // extract our new md tag and cigar
    val (cigar, md) = ObservationOperator.makeCigarAndMD(collapsedOperators)

    AlignmentRecord.newBuilder(read)
      .setCigar(cigar)
      .setMismatchingPositions(md)
      .build()
  }

  /**
   * Tests to see if the read is a candidate for realignment.
   *
   * @param alignment Alignment blocks to examine.
   * @return True if the read should be realigned.
   */
  private[realigner] def isRealignmentCandidate(
    alignment: Iterable[RealignmentBlock]): Boolean = {

    // must have at least one alignment block
    assert(alignment.nonEmpty)

    !alignment.forall(b => b match {
      case RealignableBlock(_, _) => false
      case _                      => true
    })
  }

  /**
   * Given the CIGAR and MD tag for a read, extracts the alignment blocks.
   *
   * @param read Read to extract alignment blocks for.
   * @return The alignment blocks extracted from the read.
   */
  private[realigner] def extractAlignmentOperators(
    read: AlignmentRecord): Iterable[ObservationOperator] = {

    assert(read.getReadMapped)

    if (Option(read.getCigar).fold(true)(c => c == "*") ||
      Option(read.getMismatchingPositions).isEmpty) {
      Iterable.empty
    } else {

      // parse out the Cigar and the MD tag
      val cigar = TextCigarCodec.decode(read.getCigar)
      val mdTag = MdTag(read.getMismatchingPositions,
        0L, // we're using pos as the index into the read
        cigar)

      // flatmap and parse out
      var idx = 0L
      cigar.getCigarElements().flatMap(elem => {
        val length = elem.getLength
        elem.getOperator match {
          case CigarOperator.M => {
            @tailrec def breakUpOp(opIdx: Int = 0,
                                   operators: List[ObservationOperator] = List.empty,
                                   refRun: Int = 0,
                                   ref: StringBuilder = new StringBuilder()): Iterable[ObservationOperator] = {

              def makeMismatch: List[ObservationOperator] = {
                Match(ref.length, Some(ref.toString)) :: operators
              }

              def makeMatch: List[ObservationOperator] = {
                Match(refRun) :: operators
              }

              // have we made it to the end of the operator?
              if (opIdx == length) {

                // we either end in a match or mismatch
                val finalOperators = if (refRun > 0) {
                  assert(ref.isEmpty)
                  makeMatch
                } else {
                  assert(ref.nonEmpty)
                  makeMismatch
                }

                finalOperators.reverse
                  .toIterable
              } else {
                // what have we here? specifically, is this base a match or mismatch?
                val base = mdTag.mismatchedBase(idx + opIdx)

                // update for next step
                val (newRun, newRef, newOperators) = base.fold({

                  // are we at the end of a run of sequence mismatches?
                  if (ref.nonEmpty) {
                    (1, new StringBuilder(), makeMismatch)
                  } else {
                    (refRun + 1, ref, operators)
                  }
                })(b => {

                  // are we at the end of a run of sequence matches?
                  val ops = if (refRun > 0) {
                    assert(ref.isEmpty)
                    makeMatch
                  } else {
                    operators
                  }

                  (0, ref.append(b), ops)
                })

                // recurse
                breakUpOp(opIdx + 1, newOperators, newRun, newRef)
              }
            }

            // break up operators before adding to the index
            val matchOperators = breakUpOp()
            idx += length

            matchOperators
          }
          case CigarOperator.EQ => {
            idx += length
            Iterable(Match(length))
          }
          case CigarOperator.X => {
            val newOp = Iterable(Match(length,
              Some((idx until idx + length).flatMap(i => {
                val optBase = mdTag.mismatchedBase(i)
                require(optBase.isDefined,
                  "Invalid CIGAR/MD tag combo (%s/%s). Index %d is X, but has no mismatching base in MD tag.".format(
                    read.getCigar, read.getMismatchingPositions, i))
                optBase
              }).mkString)))
            idx += length
            newOp
          }
          case CigarOperator.I => {
            Iterable(Insertion(length))
          }
          case CigarOperator.D => {
            val newOp = Iterable(Deletion((idx until idx + length).flatMap(i => {
              val optBase = mdTag.deletedBase(i)
              require(optBase.isDefined,
                "Invalid CIGAR/MD tag combo (%s/%s). Index %d is D, but has no deleted base in MD tag.".format(
                  read.getCigar, read.getMismatchingPositions, i))
              optBase
            }).mkString))
            idx += length
            newOp
          }
          case CigarOperator.S => {
            Iterable(Clipped(length))
          }
          case CigarOperator.H => {
            Iterable(Clipped(length, soft = false))
          }
          case _ => {
            throw new IllegalArgumentException("Unsupported CIGAR element %s in tag %s.".format(elem,
              read.getCigar))
          }
        }
      })
    }
  }
}
