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

import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.models.{
  Clipped,
  Deletion,
  Insertion,
  Match,
  ObservationOperator
}
import scala.annotation.tailrec
import scala.collection.mutable.StringBuilder

/**
 * Singleton object for creating collections of realignment blocks.
 */
private[realigner] object RealignmentBlock {

  /**
   * Backtracks into a previous canonical block and extracts operators.
   *
   * When we encounter a block that needs to be realigned, we need to go into
   * the previous block to extract the adjacent flanking canonical function.
   * This function does that. HERE BE DRAGONS.
   *
   * @param iter Iterator across the alignment blocks in the previous block.
   * @param kmerLength The length of flanking sequence we need.
   * @param bases The bases from the last block.
   * @param baseBuffer The bases we are building into this block. This should
   *   be empty when called unless the block we are building from is starts
   *   with an insertion. This should be reversed.
   * @param blocks The current realignment blocks.
   * @param ops The alignment operators we are building into the new block.
   * @return Returns a tuple with (the bases in the new block, the alignment
   *   operators in this block, the realignment blocks seen).
   */
  @tailrec private def backtrackBlocks(iter: Iterator[ObservationOperator],
                                       kmerLength: Int,
                                       bases: StringBuilder,
                                       baseBuffer: StringBuilder,
                                       blocks: List[RealignmentBlock],
                                       ops: List[ObservationOperator] = List.empty): (StringBuilder, List[ObservationOperator], List[RealignmentBlock]) = {

    // TODO: ops should be a ListBuffer

    // are we out of operators?
    if (!iter.hasNext) {
      (baseBuffer.reverseContents(),
        ops.reverse,
        blocks)
    } else {
      // pop from the head of the list
      val op = iter.next

      val optReturnValue: Option[(StringBuilder, List[ObservationOperator], List[RealignmentBlock])] = op match {
        case Match(length, optBases) => {
          if (length > kmerLength && optBases.isEmpty) {
            Some((baseBuffer.append(bases.takeRight(kmerLength).reverse).reverseContents(),
              (Match(kmerLength) :: ops).reverse,
              CanonicalBlock((Match(length - kmerLength) :: iter.toList).reverse) :: blocks))
          } else if (length == kmerLength && optBases.isEmpty) {
            Some((baseBuffer.append(bases.takeRight(kmerLength).reverse).reverseContents(),
              (op :: ops).reverse,
              if (iter.hasNext) {
                CanonicalBlock(iter.toList.reverse) :: blocks
              } else {
                blocks
              }))
          } else {
            None
          }
        }
        case _ => throw new IllegalStateException("Shouldn't see a non-match operator: " + op)
      }

      // alas, one cannot do optReturnValue.getOrElse here and stay tailrec
      if (optReturnValue.isDefined) {
        optReturnValue.get
      } else {
        backtrackBlocks(iter,
          kmerLength,
          bases.dropRight(op.size),
          baseBuffer.append(bases.takeRight(op.size).reverse),
          blocks,
          op :: ops)
      }
    }
  }

  /**
   * Segments an input read into canonical and realignable blocks.
   *
   * A read can contain clipped, canonical, and realignable blocks. By
   * inspection, we can determine if the ends of the read are clipped. A
   * canonical block can be determined through inspection; if we have a
   * block of matches (possibly containing simple SNP/MNPs), we label that
   * as canonical. If the block is not canonical, we label it realignable.
   *
   * @param read Read to examine.
   * @param alignment The alignment blocks for this read.
   * @param kmerLength The length k of the k-mers.
   * @return Returns a collection of realignment blocks.
   */
  def apply(read: String,
            alignment: Iterable[ObservationOperator],
            kmerLength: Int): Iterable[RealignmentBlock] = ExtractingRealignmentBlocks.time {

    @tailrec def split(iter: Iterator[ObservationOperator],
                       readBases: String,
                       baseBuffer: StringBuilder = new StringBuilder(),
                       lastBlockWasRealignable: Boolean = false,
                       lastBlockWasClipped: Boolean = false,
                       opBuffer: List[ObservationOperator] = List.empty,
                       blocks: List[RealignmentBlock] = List.empty): Iterable[RealignmentBlock] = {

      def processClip: RealignmentBlock = {
        // if we just clipped, we should only have one operator in the last
        // block, and no buffered read bases
        assert(opBuffer.size == 1)
        assert(baseBuffer.isEmpty)
        opBuffer.head match {
          case Clipped(size, soft) => ClippedBlock(size, soft = soft)
          case _ => {
            throw new IllegalStateException("Got unclipped block in clipped section.")
          }
        }
      }

      // have we run out of blocks?
      if (!iter.hasNext) {

        // we should be out of bases
        assert(readBases.isEmpty)

        // finalize the last block
        val finalBlocks = if (opBuffer.nonEmpty) {
          val finalBlock = if (lastBlockWasRealignable) {
            assert(!lastBlockWasClipped)
            RealignableBlock(baseBuffer.toString, opBuffer.reverse)
          } else if (lastBlockWasClipped) {
            processClip
          } else {
            CanonicalBlock(opBuffer.reverse)
          }
          (finalBlock :: blocks)
        } else {
          blocks
        }

        // and, return
        finalBlocks.reverse.toIterable
      } else {

        // pop from head 
        val lastOp = iter.next

        // we commonly pull the read bases off
        def newReadBases: String = readBases.drop(lastOp.size)

        def startRealignableBlock(op: ObservationOperator,
                                  readPostIndel: String,
                                  blockBases: StringBuilder): (String, StringBuilder, Boolean, Boolean, List[ObservationOperator], List[RealignmentBlock]) = {

          // call to tail recursive function that backtracks into the prior block
          val (realignableBases,
            realignableOps,
            newBlocks) = backtrackBlocks(opBuffer.toIterator,
            kmerLength,
            baseBuffer,
            blockBases,
            blocks)

          (readPostIndel,
            realignableBases,
            true,
            false,
            op :: realignableOps,
            newBlocks)
        }

        // and test
        val (newBases,
          newBaseBuffer,
          isRealignable,
          isClipped,
          newOpBuffer,
          newBlocks): (String, StringBuilder, Boolean, Boolean, List[ObservationOperator], List[RealignmentBlock]) = lastOp match {
          case Match(length, optRef) => {
            if (lastBlockWasRealignable) {

              // if we are currently in a realignable block, we can close
              // out the realignable block iff:
              //
              // * we are a sequence match
              // * we are longer than the alignment k-mer length
              if (optRef.isEmpty && length >= kmerLength) {

                // the first kmerLength bases of the match go into the
                // realignable block
                val newOpLength = length - kmerLength

                // if we have a block that is the exact length of the flank,
                // we don't have to split the last operator, it goes fully
                // in the realignable block
                val (newBB, newOB) = if (newOpLength > 0) {
                  (new StringBuilder(readBases.drop(kmerLength)
                    .take(newOpLength)),
                    List(Match(newOpLength)))
                } else {
                  (new StringBuilder(), List.empty[ObservationOperator])
                }

                (newReadBases,
                  newBB,
                  false,
                  false,
                  newOB,
                  RealignableBlock(baseBuffer.append(readBases.take(kmerLength))
                    .toString,
                    (Match(kmerLength) :: opBuffer).reverse) :: blocks)
              } else {
                (newReadBases,
                  baseBuffer.append(readBases.take(length)),
                  true,
                  false,
                  lastOp :: opBuffer,
                  blocks)
              }
            } else {

              // if the last operator was a clip, clean that up, prepend it
              // to the block list, and clean up the op buffer
              val (blks, buf) = if (lastBlockWasClipped) {
                // clipping should only occur at the ends of a read, thus
                // when the block buffer is empty or contains a hard clip block
                require(blocks.forall(blk => blk match {
                  case ClippedBlock(_, false) => true
                  case _                      => false
                }),
                  "Saw clip operator in middle of a read for operators %s".format(
                    alignment.mkString))
                (processClip :: blocks,
                  List.empty[ObservationOperator])
              } else {
                (blocks, opBuffer)
              }

              (newReadBases,
                baseBuffer.append(readBases.take(length)),
                false,
                false,
                lastOp :: buf,
                blks)
            }
          }
          case Clipped(length, soft) => {
            // hard clipped bases aren't in the read, so don't drop bases
            // from the read if we are hard clipped
            def dropIfSoft: String = {
              if (soft) {
                readBases.drop(length)
              } else {
                readBases
              }
            }

            // when the last block was clipped, we need to process/validate
            // the clip, and we can reuse the empty base buffer
            val (buffer, blks) = if (lastBlockWasClipped) {
              // two consecutive clip blocks are legal if:
              //
              // * we are at the start of the read and the first block is
              //   hard clipped and the second is soft clipped
              // * we are at the end of the read and the last block is hard
              //   clipped and the previous block is soft clipped
              //
              // we can't validate that we are at the end of the read
              // here---we do that in the Match case---but we can validate
              // all of the other positional requirements
              assert(opBuffer.nonEmpty)
              opBuffer.head match {
                case Clipped(_, isSoft) => {
                  require((!isSoft && soft && blocks.isEmpty) ||
                    (isSoft && !soft),
                    "Illegal sequence of clipped bases in %s.".format(alignment.mkString))
                }
                case _ => assert(false) // should not see non-clipped operator
              }

              // also, the base buffer should be empty after a clip
              assert(baseBuffer.isEmpty)

              (baseBuffer,
                processClip :: blocks)
            } else if (lastBlockWasRealignable) {
              (new StringBuilder(),
                RealignableBlock(baseBuffer.toString, opBuffer.reverse) :: blocks)
            } else if (opBuffer.nonEmpty) {
              (new StringBuilder(),
                CanonicalBlock(opBuffer.reverse) :: blocks)
            } else {
              (new StringBuilder(),
                blocks)
            }

            (dropIfSoft,
              buffer,
              false,
              true,
              List(lastOp),
              blks)
          }
          case Insertion(length) => {
            if (lastBlockWasClipped) {
              (newReadBases,
                new StringBuilder(readBases.take(length)),
                true,
                false,
                List(lastOp),
                processClip :: blocks)
            } else if (!lastBlockWasRealignable) {
              startRealignableBlock(lastOp,
                readBases.drop(length),
                new StringBuilder(readBases.take(length)).reverse)
            } else {
              (newReadBases,
                baseBuffer.append(readBases.take(length)),
                true,
                false,
                lastOp :: opBuffer,
                blocks)
            }
          }
          case Deletion(_) => {
            if (lastBlockWasClipped) {
              (readBases,
                new StringBuilder(),
                true,
                false,
                List(lastOp),
                processClip :: blocks)
            } else if (!lastBlockWasRealignable) {
              startRealignableBlock(lastOp,
                readBases,
                new StringBuilder())
            } else {
              (readBases,
                baseBuffer,
                true,
                false,
                lastOp :: opBuffer,
                blocks)
            }
          }
        }

        split(iter,
          newBases,
          newBaseBuffer,
          isRealignable,
          isClipped,
          newOpBuffer,
          newBlocks)
      }
    }

    split(alignment.toIterator,
      read)
  }
}

/**
 * Trait for describing a block possibly for realignment.
 */
private[realigner] sealed trait RealignmentBlock {

  /**
   * Folds over this block and optionally applies a realignment function.
   *
   * @param fn Realignment function to optionally apply.
   * @return Returns alignment blocks. The realignment function is applied if
   *   this block is not clipped nor canonical.
   */
  def fold(fn: (String, Iterable[ObservationOperator]) => Iterable[ObservationOperator]): Iterable[ObservationOperator]
}

/**
 * Represents a block that was clipped.
 *
 * @param length The number of bases that were clipped.
 * @param soft If true, the bases were soft clipped.
 */
private[realigner] case class ClippedBlock(
    length: Int, soft: Boolean = true) extends RealignmentBlock {

  /**
   * @param fn Realignment function to optionally apply.
   * @return Returns a clipped alignment block. The realignment function is not
   *   applied.
   */
  def fold(fn: (String, Iterable[ObservationOperator]) => Iterable[ObservationOperator]): Iterable[ObservationOperator] = {
    Iterable(Clipped(length, soft = soft))
  }
}

/**
 * Represents a block that was canonically aligned.
 *
 * @param alignments The alignment operators for this block.
 */
private[realigner] case class CanonicalBlock(
    alignments: Iterable[ObservationOperator]) extends RealignmentBlock {

  /**
   * @param fn Realignment function to optionally apply.
   * @return Returns the alignment operators for this block. The realignment
   *   function is not applied.
   */
  def fold(fn: (String, Iterable[ObservationOperator]) => Iterable[ObservationOperator]): Iterable[ObservationOperator] = {
    alignments
  }
}

/**
 * Represents a block where the original alignment may not be canonical.
 *
 * @param read The read sequence to realign.
 * @param alignments The present alignment of this block.
 */
private[realigner] case class RealignableBlock(read: String,
                                               alignments: Iterable[ObservationOperator]) extends RealignmentBlock {

  /**
   * @param fn Realignment function to optionally apply.
   * @return Returns the alignment operators for this block. The realignment
   *   function is applied to the read string and the alignments for this block.
   */
  def fold(fn: (String, Iterable[ObservationOperator]) => Iterable[ObservationOperator]): Iterable[ObservationOperator] = {
    fn(read, alignments)
  }
}
