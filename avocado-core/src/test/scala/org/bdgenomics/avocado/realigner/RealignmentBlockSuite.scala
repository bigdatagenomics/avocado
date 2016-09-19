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
  Clipped,
  Deletion,
  Match,
  Insertion,
  ObservationOperator
}
import org.scalatest.FunSuite

object RealignmentBlockSuite {

  def foldFn(bases: String,
             ops: Iterable[ObservationOperator]): Iterable[ObservationOperator] = {
    assert(bases == "ACGT")
    assert(ops.nonEmpty)
    Iterable(Match(bases.length, Some(bases)))
  }
}

class RealignmentBlockSuite extends FunSuite {

  test("folding over a clip returns the clip operator, soft clip") {
    val ops = ClippedBlock(5).fold(RealignmentBlockSuite.foldFn)
      .toSeq

    assert(ops.size === 1)
    assert(ops(0) === Clipped(5))
  }

  test("folding over a clip returns the clip operator, hard clip") {
    val ops = ClippedBlock(5, false).fold(RealignmentBlockSuite.foldFn)
      .toSeq

    assert(ops.size === 1)
    assert(ops(0) === Clipped(5, false))
  }

  test("folding over a canonical block returns the original alignment") {
    val ops = CanonicalBlock(Iterable(Match(3),
      Match(1, Some("A")),
      Match(4)))
      .fold(RealignmentBlockSuite.foldFn)
      .toSeq

    assert(ops.size === 3)
    assert(ops(0) === Match(3))
    assert(ops(1) === Match(1, Some("A")))
    assert(ops(2) === Match(4))
  }

  test("violate an invariant of the fold function, part 1") {
    intercept[AssertionError] {
      RealignableBlock("TGCA", Iterable(Match(4))).fold(RealignmentBlockSuite.foldFn)
    }
  }

  test("violate an invariant of the fold function, part 2") {
    intercept[AssertionError] {
      RealignableBlock("ACGT", Iterable()).fold(RealignmentBlockSuite.foldFn)
    }
  }

  test("apply the fold function on a realignable block") {
    val ops = RealignableBlock("ACGT", Iterable(Match(4)))
      .fold(RealignmentBlockSuite.foldFn)
      .toSeq

    assert(ops.size === 1)
    assert(ops(0) === Match(4, Some("ACGT")))
  }

  test("having a clip in the middle of a read is illegal") {
    intercept[IllegalArgumentException] {
      RealignmentBlock("TTTAAGGG",
        Iterable(Match(3), Clipped(2), Match(3)),
        3)
    }
  }

  test("can't have two soft clips back to back") {
    intercept[IllegalArgumentException] {
      RealignmentBlock("TTTAAGGG",
        Iterable(Clipped(3), Clipped(2), Match(3)),
        3)
    }
  }

  test("a read that is an exact sequence match is canonical") {
    val blocks = RealignmentBlock("ACTGTTGTTACAGA",
      Iterable(Match(14)),
      5).toSeq

    assert(blocks.size === 1)
    blocks(0) match {
      case CanonicalBlock(alignments) => {
        val aln = alignments.toSeq
        assert(aln.size === 1)
        assert(aln(0) === Match(14))
      }
      case _ => assert(false)
    }
  }

  def zipAndCheckAlignments(alignments: Iterable[ObservationOperator],
                            ops: Iterable[ObservationOperator]) = {
    assert(alignments.size === ops.size)
    alignments.zip(ops).foreach(p => assert(p._1 === p._2))
  }

  test("hard clip before soft clip is ok at start of read") {
    val blocks = RealignmentBlock("TAAGG",
      Iterable(Clipped(3, soft = false), Clipped(2), Match(3)),
      3).toSeq

    assert(blocks.size === 3)
    assert(blocks(0) === ClippedBlock(3, soft = false))
    assert(blocks(1) === ClippedBlock(2))
    blocks(2) match {
      case CanonicalBlock(alignments) => {
        assert(alignments.size === 1)
        assert(alignments.head === Match(3))
      }
      case _ => assert(false)
    }
  }

  test("hard clip after soft clip is ok at end of read") {
    val blocks = RealignmentBlock("TTTAAG",
      Iterable(Match(3), Clipped(3), Clipped(2, soft = false)),
      3).toSeq

    assert(blocks.size === 3)
    blocks(0) match {
      case CanonicalBlock(alignments) => {
        assert(alignments.size === 1)
        assert(alignments.head === Match(3))
      }
      case _ => assert(false)
    }
    assert(blocks(1) === ClippedBlock(3))
    assert(blocks(2) === ClippedBlock(2, soft = false))
  }

  test("a read with a single snp is canonical") {
    val ops = Iterable(Match(6),
      Match(1, Some("G")),
      Match(7))
    val blocks = RealignmentBlock("ACTGTTATTACAGA",
      ops,
      5).toSeq

    assert(blocks.size === 1)
    blocks(0) match {
      case CanonicalBlock(alignments) => {
        zipAndCheckAlignments(alignments, ops)
      }
      case _ => assert(false)
    }
  }

  test("a read containing an indel with exact flanks is wholly realignable") {
    val readSeq = "ACTGTTATTACAGA"
    val ops = Iterable(Match(5),
      Insertion(4),
      Match(5))
    val blocks = RealignmentBlock(readSeq,
      ops,
      5).toSeq

    assert(blocks.size === 1)
    blocks(0) match {
      case RealignableBlock(bases, alignments) => {
        assert(bases === readSeq)
        zipAndCheckAlignments(alignments, ops)
      }
      case _ => assert(false)
    }
  }

  test("a read containing an indel with exact flanks is wholly realignable, with soft clipped bases") {
    val readSeq = "ACTGTTATTACAGA"
    val ops1 = Iterable(Clipped(4))
    val ops2 = Iterable(Match(5),
      Deletion("AC"),
      Match(5))
    val blocks = RealignmentBlock(readSeq,
      ObservationOperator.collapse(Iterable.concat(ops1, ops2)),
      5).toSeq

    assert(blocks.size === 2)
    blocks(0) match {
      case ClippedBlock(length, true) => {
        assert(length === 4)
      }
      case _ => assert(false)
    }
    blocks(1) match {
      case RealignableBlock(bases, alignments) => {
        assert(bases === readSeq.drop(4))
        zipAndCheckAlignments(alignments, ops2)
      }
      case _ => assert(false)
    }
  }

  test("a read containing an indel with longer flanks can be split into multiple blocks") {
    val readSeq = "TTACACTGTTATTACATG"
    val ops1 = Iterable(Match(4))
    val ops2 = Iterable(Match(5),
      Insertion(2),
      Match(5))
    val ops3 = Iterable(Clipped(2))
    val blocks = RealignmentBlock(readSeq,
      ObservationOperator.collapse(Iterable.concat(ops1, ops2, ops3)),
      5).toSeq

    assert(blocks.size === 3)
    blocks(0) match {
      case CanonicalBlock(alignments) => {
        zipAndCheckAlignments(alignments, ops1)
      }
      case _ => assert(false)
    }
    blocks(1) match {
      case RealignableBlock(bases, alignments) => {
        assert(bases === readSeq.drop(4).dropRight(2))
        zipAndCheckAlignments(alignments, ops2)
      }
      case _ => assert(false)
    }
    blocks(2) match {
      case ClippedBlock(length, true) => {
        assert(length === 2)
      }
      case _ => assert(false)
    }
  }

  test("a read containing an indel with longer flanks on both sides can be split into multiple blocks") {
    val readSeq = "TTACACTGTTATTACATGG"
    val ops1 = Iterable(Clipped(10, soft = false))
    val ops2 = Iterable(Match(4))
    val ops3 = Iterable(Match(5),
      Insertion(2),
      Match(5))
    val ops4 = Iterable(Match(1, Some("A")), Match(2))
    val blocks = RealignmentBlock(readSeq,
      ObservationOperator.collapse(Iterable.concat(ops1, ops2, ops3, ops4)),
      5).toSeq

    assert(blocks.size === 4)
    blocks(0) match {
      case ClippedBlock(length, false) => {
        assert(length === 10)
      }
      case _ => assert(false)
    }
    blocks(1) match {
      case CanonicalBlock(alignments) => {
        zipAndCheckAlignments(alignments, ops2)
      }
      case _ => assert(false)
    }
    blocks(2) match {
      case RealignableBlock(bases, alignments) => {
        assert(bases === readSeq.drop(4).dropRight(3))
        zipAndCheckAlignments(alignments, ops3)
      }
      case _ => assert(false)
    }
    blocks(3) match {
      case CanonicalBlock(alignments) => {
        zipAndCheckAlignments(alignments, ops4)
      }
      case _ => assert(false)
    }
  }
}
