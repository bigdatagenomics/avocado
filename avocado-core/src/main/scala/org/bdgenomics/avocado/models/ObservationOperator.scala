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
package org.bdgenomics.avocado.models

import htsjdk.samtools.{ CigarOperator, TextCigarCodec }
import org.bdgenomics.adam.util.MdTag
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.formats.avro.AlignmentRecord
import scala.annotation.tailrec
import scala.collection.JavaConversions._
import scala.collection.mutable.StringBuilder

/**
 * Companion object for ObservationOperator trait.
 *
 * Implements convenience operations on top of collections of
 * ObservationOperators.
 */
object ObservationOperator {

  /**
   * Given the CIGAR and MD tag for a read, extracts the alignment blocks.
   *
   * @param read Read to extract alignment blocks for.
   * @return The alignment blocks extracted from the read.
   */
  private[avocado] def extractAlignmentOperators(
    read: AlignmentRecord): Iterable[ObservationOperator] = ExtractingAlignment.time {

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

  /**
   * Collapses down a collection of ObservationOperators.
   *
   * Where two blocks that are the same type of operator are seen, merges the
   * blocks together. Also, filters out any 0 length blocks.
   *
   * @param ops The collection of operators to collapse down.
   * @return A new collection of operators where adjacent identical operators
   *   have been merged together.
   */
  private[avocado] def collapse(ops: Iterable[ObservationOperator]): Iterable[ObservationOperator] = {

    @tailrec def collapser(
      iter: Iterator[ObservationOperator],
      head: Option[ObservationOperator] = None,
      collapsed: List[ObservationOperator] = List.empty): Iterable[ObservationOperator] = {
      if (!iter.hasNext) {
        val finalList = if (head.isDefined) {
          head.get :: collapsed
        } else {
          collapsed
        }

        finalList.reverse.toIterable
      } else {
        val next = iter.next
        val (newHead, newList): (ObservationOperator, List[ObservationOperator]) = head.fold((next, collapsed))(h => {
          (next, h) match {
            case (a: Match, b: Match) => {
              if (a.isSequenceMatch && b.isSequenceMatch) {
                (Match(a.size + b.size), collapsed)
              } else if (a.isSequenceMismatch && b.isSequenceMismatch) {
                // note: b needs to go ahead of a
                (Match(b.size + a.size,
                  Some(b.optMismatchingBases.get + a.optMismatchingBases.get)),
                  collapsed)
              } else {
                (a, b :: collapsed)
              }
            }
            case (Insertion(a), Insertion(b))     => (Insertion(a + b), collapsed)
            case (Deletion(aDel), Deletion(bDel)) => (Deletion(bDel + aDel), collapsed)
            case _                                => (next, h :: collapsed)
          }
        })

        collapser(iter, head = Some(newHead), collapsed = newList)
      }
    }

    collapser(ops.filter(_.nonEmpty).toIterator)
  }

  /**
   * Extracts a reference sequence from a read and alignment operators.
   *
   * @param read Read sequence to extract reference from.
   * @param ops The alignment operators to use for extracting the reference.
   * @return The reference sequence covering these blocks.
   */
  private[avocado] def extractReference(
    read: String,
    ops: Iterable[ObservationOperator]): String = ExtractingReference.time {

    def errMsg: String = {
      "%s, %s".format(ops.mkString(","), read)
    }

    @tailrec def extract(iter: Iterator[ObservationOperator],
                         trimmedRead: String,
                         reference: StringBuilder = new StringBuilder()): String = {

      if (!iter.hasNext) {
        require(trimmedRead.isEmpty,
          "Processed alignment operators and read (%s) but still have bases left: %s".format(
            errMsg, trimmedRead))
        reference.toString
      } else {

        // pop from the iterator
        val op = iter.next

        def trimOnly(l: Int): String = {
          require(trimmedRead.length > l,
            "Got %s but don't have enough bases (%s) to drop for read (%s).".format(
              op, trimmedRead, errMsg))
          trimmedRead.drop(l)
        }

        // and test
        val newRead = op match {
          case Clipped(length, true) => trimOnly(length)
          case Insertion(length)     => trimOnly(length)
          case Clipped(_, false) => {
            // hard clip is a no-op
            trimmedRead
          }
          case Deletion(ref) => {
            // append deleted sequence
            reference.append(ref)

            // don't do anything to read
            trimmedRead
          }
          case Match(length, optRef) => {
            // append sequence to ref
            reference.append(optRef.getOrElse(trimmedRead.take(length)))

            // trim front of read
            trimmedRead.drop(length)
          }
        }

        // recurse
        extract(iter, newRead, reference)
      }
    }

    extract(ops.toIterator, read)
  }

  /**
   * Turns a collection of operators into a text cigar and MD tag.
   *
   * @param ops Alignment operators to turn into a cigar/MD tag.
   * @return Returns a tuple of (CIGAR, MD tag).
   */
  private[avocado] def makeCigarAndMD(ops: Iterable[ObservationOperator]): (String, String) = {

    @tailrec def make(iter: Iterator[ObservationOperator],
                      matchRun: Int = 0,
                      cigar: StringBuilder = new StringBuilder(),
                      md: StringBuilder = new StringBuilder()): (String, String) = {
      if (!iter.hasNext) {

        // add the final run to the md tag
        md.append(matchRun)

        // build our string builders
        (cigar.toString, md.toString)
      } else {

        // pop from the head of the iterator and append to the cigar
        val op = iter.next
        cigar.append(op.toString)

        // test and modify the cigar
        val newMatchRun = op match {
          case Match(length, None) => {
            matchRun + length
          }
          case Match(_, Some(ref)) => {

            // append the number of bases that just matched
            md.append(matchRun)

            // append the first ref base
            md.append(ref.head)

            // loop over the remaining ref bases and add
            ref.tail.foreach(base => {
              md.append(0)
              md.append(base)
            })

            // a mismatch resets the match length
            0
          }
          case Deletion(ref) => {

            // append the number of bases that just matched
            md.append(matchRun)

            // append the deletion caret
            md.append('^')

            // append the deleted sequence
            md.append(ref)

            // a deletion resets the match length
            0
          }
          case _ => {

            // we are in an insertion or a clip
            matchRun
          }
        }

        make(iter, newMatchRun, cigar, md)
      }
    }

    make(ops.toIterator)
  }
}

/**
 * A trait that implements an alignment operator.
 *
 * These operators represent a limited subset of the CIGAR operators:
 *
 * * Alignment Match (M --> X/=)
 * * Insertion (I)
 * * Deletion (D)
 * * Clip (S/H)
 */
sealed trait ObservationOperator {

  /**
   * The length of this alignment block.
   */
  val size: Int

  /**
   * @return True if this block has non-zero length.
   */
  def nonEmpty: Boolean = size > 0
}

/**
 * An alignment match block (M), further broken down into a sequence (mis)match.
 *
 * @param size The size of this block.
 * @param optMismatchingBases If this block is a sequence match (=), this will
 *   be empty. If this block is a sequence mismatch (X), this will contain the
 *   mismatching reference bases.
 */
case class Match(size: Int,
                 optMismatchingBases: Option[String] = None) extends ObservationOperator {

  optMismatchingBases.foreach(r => require(r.length == size,
    "Ref bases must have same length as block."))

  /**
   * @return True if there are no mismatching bases.
   */
  def isSequenceMatch: Boolean = optMismatchingBases.isEmpty

  /**
   * @return True if there are mismatching bases.
   */
  def isSequenceMismatch: Boolean = optMismatchingBases.isDefined

  /**
   * @return Returns "<length>=" if a match, "<length>X" if a mismatch.
   */
  override def toString: String = {
    if (optMismatchingBases.isEmpty) {
      "%d=".format(size)
    } else {
      "%dX".format(size)
    }
  }
}

/**
 * A block representing an insertion into the reference.
 *
 * @param size The insertion length.
 */
case class Insertion(size: Int) extends ObservationOperator {

  /**
   * @return Returns "<length>I".
   */
  override def toString: String = "%dI".format(size)
}

/**
 * A block representing a deletion from the reference.
 *
 * @param basesDeleted The bases deleted from the reference.
 */
case class Deletion(basesDeleted: String) extends ObservationOperator {

  /**
   * The number of bases deleted from the reference.
   */
  val size = basesDeleted.size

  /**
   * @return Returns "<length>D".
   */
  override def toString: String = "%dD".format(size)
}

/**
 * A block representing bases clipped from the read.
 *
 * @param size The number of bases clipped.
 * @param soft True if soft clipped, false if hard clipped.
 */
case class Clipped(size: Int,
                   soft: Boolean = true) extends ObservationOperator {

  /**
   * @return Returns "<length>S" if soft clipped, "<length>H" if hard clipped.
   */
  override def toString: String = {
    "%d%s".format(size, if (soft) "S" else "H")
  }
}
