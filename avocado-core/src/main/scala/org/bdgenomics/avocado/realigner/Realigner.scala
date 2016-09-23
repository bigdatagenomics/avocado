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

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.avocado.Timers._
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
        ProcessingReadForRealignment.time {

          // convert the cigar and MD tag into alignment operators
          val alignment = ObservationOperator.extractAlignmentOperators(r)

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
              try {
                realignRead(r, realignmentBlocks, kmerLength)
              } catch {
                case t: Throwable => {
                  log.warn("Realigning %s failed with exception %s.".format(
                    r.getReadName, t))
                  r
                }
              }
            } else {
              r
            }
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
    kmerLength: Int): AlignmentRecord = RealigningRead.time {

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

    FinalizingRealignment.time {
      // we may have a run of likewise blocks, so collapse on down
      val collapsedOperators = ObservationOperator.collapse(realignedOperators)

      // extract our new md tag and cigar
      val (cigar, md) = ObservationOperator.makeCigarAndMD(collapsedOperators)

      AlignmentRecord.newBuilder(read)
        .setCigar(cigar)
        .setMismatchingPositions(md)
        .build()
    }
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
}
