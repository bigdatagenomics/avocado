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

import org.bdgenomics.adam.models.{
  SequenceDictionary,
  SequenceRecord,
  RecordGroup,
  RecordGroupDictionary
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.avocado.models.{
  Clipped,
  Deletion,
  Insertion,
  Match,
  ObservationOperator
}
import org.bdgenomics.formats.avro.AlignmentRecord

class RealignerSuite extends AvocadoFunSuite {

  test("realignment candidate code needs at least one block") {
    intercept[AssertionError] {
      Realigner.isRealignmentCandidate(Iterable.empty)
    }
  }

  test("read is not a realignment candidate if it is canonical") {
    assert(!Realigner.isRealignmentCandidate(Iterable(
      CanonicalBlock(Iterable(Match(20))))))
  }

  test("read is not a realignment candidate if it is canonical and clipped") {
    assert(!Realigner.isRealignmentCandidate(Iterable(
      CanonicalBlock(Iterable(Match(20),
        Clipped(10))))))
  }

  test("read is a realignment candidate if there is at least one non-canonical block") {
    assert(Realigner.isRealignmentCandidate(Iterable(
      RealignableBlock("ACATT", Iterable(
        Match(2),
        Insertion(1),
        Match(2))))))
  }

  def realign(read: AlignmentRecord,
              kmerLength: Int): AlignmentRecord = {

    val alignment = ObservationOperator.extractAlignmentOperators(read)

    val ops = RealignmentBlock(read.getSequence,
      alignment,
      kmerLength)

    Realigner.realignRead(read, ops, kmerLength)
  }

  test("realign an indel that is not left normalized") {
    val read = AlignmentRecord.newBuilder()
      .setReadName("A_READ")
      .setReadMapped(true)
      .setStart(10L)
      .setEnd(21L)
      .setSequence("ATACACTTAG")
      .setCigar("4M2I4M")
      .setMismatchingPositions("8")
      .build()

    val newAlignment = realign(read, 4)

    assert(newAlignment.getCigar === "2=2I6=")
    assert(newAlignment.getMismatchingPositions === "8")
  }

  test("realign a mnp expressed as a complex indel") {
    val read = AlignmentRecord.newBuilder()
      .setReadName("A_READ")
      .setReadMapped(true)
      .setStart(10L)
      .setEnd(33L)
      .setSequence("ACTCTGAAACTCCAACCACTGGA")
      .setCigar("10M3I3D10M")
      .setMismatchingPositions("10^AGA10")
      .build()

    val newAlignment = realign(read, 10)

    assert(newAlignment.getCigar === "10=3X10=")
    assert(newAlignment.getMismatchingPositions === "10A0G0A10")
  }

  test("realign two snps expressed as a complex indel") {
    val read = AlignmentRecord.newBuilder()
      .setReadName("A_READ")
      .setReadMapped(true)
      .setStart(10L)
      .setEnd(34L)
      .setSequence("ACTCTGAAACTGCAACCACTGGA")
      .setCigar("10M3I3D10M")
      .setMismatchingPositions("10^AGA10")
      .build()

    val newAlignment = realign(read, 10)

    assert(newAlignment.getCigar === "10=1X1=1X10=")
    assert(newAlignment.getMismatchingPositions === "10A1A10")
  }

  test("align sequence with a complex deletion") {
    val obs = Aligner.align("ACTCTGAATAGGGAACCACTGGA",
      "ACTCTGAATAC__AACCACTGGA".filter(_ != '_'),
      10).toSeq

    assert(obs.size === 4)
    assert(obs(0) === Match(10))
    assert(obs(1) === Deletion("GG"))
    assert(obs(2) === Match(1, Some("G")))
    assert(obs(3) === Match(10))
  }

  test("realign a read with a complex deletion") {
    // alt: ACTCTGAATAC__AACCACTGGA
    // ref: ACTCTGAATAGGGAACCACTGGA
    //                |
    val read = AlignmentRecord.newBuilder()
      .setReadName("A_READ")
      .setReadMapped(true)
      .setStart(10L)
      .setEnd(34L)
      .setSequence("ACTCTGAATACAACCACTGGA")
      .setCigar("11M2D10M")
      .setMismatchingPositions("10G0^GG10")
      .build()

    val newAlignment = realign(read, 10)

    assert(newAlignment.getCigar === "10=2D1X10=")
    assert(newAlignment.getMismatchingPositions === "10^GG0G10")
  }

  test("realign a read with a snp and deletion separated by a flank") {
    val read = AlignmentRecord.newBuilder()
      .setReadName("A_READ")
      .setReadMapped(true)
      .setStart(10L)
      .setEnd(43L)
      .setSequence("ACTCTGAAACTGCCTAGCACCAAACCACTGGA")
      .setCigar("24M1D8M")
      .setMismatchingPositions("10A13^A8")
      .build()

    val newAlignment = realign(read, 10)

    assert(newAlignment.getCigar === "10=1X10=1D11=")
    assert(newAlignment.getMismatchingPositions === "10A10^A11")
  }

  test("realigning a repetative read will fire an assert") {
    val read = AlignmentRecord.newBuilder()
      .setReadName("A_READ")
      .setReadMapped(true)
      .setStart(10L)
      .setEnd(17L)
      .setSequence("TCAAAAAAGG")
      .setCigar("3M4I3M")
      .setMismatchingPositions("6")
      .build()

    intercept[AssertionError] {
      realign(read, 3)
    }
  }

  def makeAndRealignRdd(reads: Seq[AlignmentRecord],
                        kmerLength: Int): Array[AlignmentRecord] = {
    val gRdd = AlignmentRecordRDD(sc.parallelize(reads),
      SequenceDictionary(SequenceRecord("ctg", 50L)),
      RecordGroupDictionary(Seq(RecordGroup("rg", "rg"))))

    // realign the genomic rdd
    val realignedRdd = Realigner.realign(gRdd, kmerLength)

    // collect the reads
    realignedRdd.rdd.collect()
  }

  sparkTest("realign a set of reads around an insert") {
    // insertion sequence:
    // ins: AATGAGACTTACATCATTAAAACCGTGTGGACACA
    // ref: AATGAGACTTACATCATTAA__CCGTGTGGACACA
    val sequence = "AATGAGACTTACATCATTAAAACCGTGTGGACACA"
    val insertStart = 20
    val readLength = insertStart + 6 + 2

    // generate 7 reads with a 6bp flank
    val reads = (0 until 7).map(rId => {
      val basesBeforeInsert = insertStart - rId
      val basesAfterInsert = 6 + rId

      AlignmentRecord.newBuilder()
        .setReadName(rId.toString)
        .setContigName("ctg")
        .setRecordGroupName("rg")
        .setReadMapped(true)
        .setSequence(sequence.drop(rId).take(readLength))
        .setStart(rId.toLong)
        .setEnd((rId + readLength - 2 + 1).toLong)
        .setCigar("%dM2I%dM".format(basesBeforeInsert, basesAfterInsert))
        .setMismatchingPositions((readLength - 2).toString)
        .build()
    })

    // make into a genomic rdd
    val newReads = makeAndRealignRdd(reads, 6)

    assert(newReads.size === 7)
    newReads.foreach(r => {
      val rId = r.getReadName.toInt

      // these values are different from above because original alignments were
      // not left justified
      val basesBeforeInsert = insertStart - rId - 2
      val basesAfterInsert = 8 + rId

      assert(r.getCigar === "%d=2I%d=".format(basesBeforeInsert, basesAfterInsert))
      assert(r.getMismatchingPositions === (readLength - 2).toString)
    })
  }

  sparkTest("realign a set of reads around a deletion") {
    // deletion sequence:
    // del: AGGTCTGAATGAGACTTA__TCATTAACCGTGTGGACACA
    // ref: AGGTCTGAATGAGACTTACATCATTAACCGTGTGGACACA
    val sequence = "AGGTCTGAATGAGACTTATCATTAACCGTGTGGACACA"
    val deleteStart = 18
    val readLength = deleteStart + 8

    // generate 10 reads with a 8bp flank
    val reads = (0 until 10).map(rId => {
      val basesBeforeDelete = deleteStart - rId
      val basesAfterDelete = 8 + rId

      AlignmentRecord.newBuilder()
        .setReadName(rId.toString)
        .setContigName("ctg")
        .setRecordGroupName("rg")
        .setReadMapped(true)
        .setSequence(sequence.drop(rId).take(readLength))
        .setStart(rId.toLong)
        .setEnd((rId + readLength + 2 + 1).toLong)
        .setCigar("%dM2D%dM".format(basesBeforeDelete, basesAfterDelete))
        .setMismatchingPositions("%d^CA%d".format(basesBeforeDelete, basesAfterDelete))
        .build()
    })

    // make into a genomic rdd
    val newReads = makeAndRealignRdd(reads, 8)

    assert(newReads.size === 10)
    newReads.foreach(r => {
      val rId = r.getReadName.toInt

      // these values are different from above because original alignments were
      // not left justified
      val basesBeforeDelete = deleteStart - rId - 1
      val basesAfterDelete = 9 + rId

      assert(r.getCigar === "%d=2D%d=".format(basesBeforeDelete, basesAfterDelete))
      assert(r.getMismatchingPositions === "%d^AC%d".format(basesBeforeDelete, basesAfterDelete))
    })
  }

  sparkTest("realigning a read with a repeat will return the original read") {
    val read = AlignmentRecord.newBuilder()
      .setReadName("A_READ")
      .setReadMapped(true)
      .setStart(10L)
      .setEnd(17L)
      .setSequence("TCAAAAAAGG")
      .setCigar("3M4I3M")
      .setMismatchingPositions("6")
      .build()

    // make into a genomic rdd
    val newReads = makeAndRealignRdd(Seq(read), 3)

    // should have one read
    assert(newReads.size === 1)
    assert(newReads.head === read)
  }

  sparkTest("one sample read should fail due to a repeat, all others should realign") {
    val readFile = resourceUrl("NA12878_reads.sam").toString
    val reads = sc.loadAlignments(readFile)
      .rdd
      .collect
      .toSeq

    var asserts = 0
    reads.foreach(r => {
      try {
        realign(r, 20)
      } catch {
        case ae: AssertionError => {
          asserts += 1
        }
      }
    })

    assert(asserts === 1)
  }
}
