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
import org.bdgenomics.adam.rdd.read.{ AlignedReadRDD, AlignmentRecordRDD }
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.AlignmentRecord

trait SparkRealignerSuite extends AvocadoFunSuite {

  def realign(rdd: AlignmentRecordRDD,
              kmerLength: Int): AlignmentRecordRDD

  def makeAndRealignRdd(reads: Seq[AlignmentRecord],
                        kmerLength: Int): Array[AlignmentRecord] = {
    val gRdd = AlignedReadRDD(sc.parallelize(reads),
      SequenceDictionary(SequenceRecord("ctg", 50L)),
      RecordGroupDictionary(Seq(RecordGroup("rg", "rg"))))

    // realign the genomic rdd
    val realignedRdd = realign(gRdd, kmerLength)

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
}
