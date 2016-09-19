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

import scala.annotation.tailrec
import scala.collection.mutable.StringBuilder

/**
 * Companion object for ObservationOperator trait.
 *
 * Implements convenience operations on top of collections of
 * ObservationOperators.
 */
object ObservationOperator {

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
  private[avocado] def extractReference(read: String,
                                        ops: Iterable[ObservationOperator]): String = {

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
