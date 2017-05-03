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
package org.bdgenomics.avocado.genotyping

import org.bdgenomics.adam.models.{
  RecordGroupDictionary,
  SequenceDictionary,
  SequenceRecord
}
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.{ AlignmentRecord, Variant }

class DiscoverVariantsSuite extends AvocadoFunSuite {

  val unalignedRead = AlignmentRecord.newBuilder()
    .setReadMapped(false)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .build

  val perfectReadMCigar = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("1")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .setCigar("8M")
    .setMismatchingPositions("8")
    .build

  val perfectReadEqCigar = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("1")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .setCigar("8=")
    .setMismatchingPositions("8")
    .build

  val snpReadMCigar = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("1")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .setCigar("8M")
    .setMismatchingPositions("4C3")
    .build

  val snpReadEqCigar = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("1")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .setCigar("4=1X3=")
    .setMismatchingPositions("4C3")
    .build

  val snpReadHardClip = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("1")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .setCigar("2H8M")
    .setMismatchingPositions("4C3")
    .build

  val snpReadSoftClip = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("1")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("TGACACATGA")
    .setQual("!!!!!!!!!!")
    .setCigar("2S8M")
    .setMismatchingPositions("4C3")
    .build

  val insertRead = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("2")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("ACACTTATGA")
    .setQual("!!!!!!!!!!")
    .setCigar("4M2I4M")
    .setMismatchingPositions("8")
    .build

  val deleteRead = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("3")
    .setStart(10L)
    .setEnd(20L)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .setCigar("4M2D4M")
    .setMismatchingPositions("4^TT4")
    .build

  val mnpRead = AlignmentRecord.newBuilder()
    .setReadMapped(true)
    .setContigName("3")
    .setStart(10L)
    .setEnd(18L)
    .setSequence("ACACATGA")
    .setQual("!!!!!!!!")
    .setCigar("8M")
    .setMismatchingPositions("3T0T3")
    .build

  test("no variants in unaligned read") {
    assert(DiscoverVariants.variantsInRead(unalignedRead, 0).isEmpty)
  }

  sparkTest("no variants in rdd with unaligned read") {
    val unalignedRdd = sc.parallelize(Seq(unalignedRead))
    assert(DiscoverVariants.variantsInRdd(unalignedRdd).count === 0)
  }

  test("no variants in read that is a perfect sequence match") {
    assert(DiscoverVariants.variantsInRead(perfectReadMCigar, 0).isEmpty)
    assert(DiscoverVariants.variantsInRead(perfectReadEqCigar, 0).isEmpty)
  }

  sparkTest("no variants in rdd with sequence match reads") {
    val matchRdd = sc.parallelize(Seq(perfectReadMCigar, perfectReadEqCigar))
    assert(DiscoverVariants.variantsInRdd(matchRdd).count === 0)
  }

  def validateSnp(snp: Variant) {
    assert(snp.getContigName() === "1")
    assert(snp.getStart() === 14L)
    assert(snp.getEnd() === 15L)
    assert(snp.getReferenceAllele === "C")
    assert(snp.getAlternateAllele === "A")
  }

  test("find snp in read with a 1bp sequence mismatch") {
    def testSnp(read: AlignmentRecord) {
      val variants = DiscoverVariants.variantsInRead(read, 0)
      assert(variants.size === 1)
      validateSnp(variants.head.toVariant)
      val highQualVariants = DiscoverVariants.variantsInRead(read, 1)
      assert(highQualVariants.isEmpty)
    }
    testSnp(snpReadMCigar)
    testSnp(snpReadEqCigar)
    testSnp(snpReadSoftClip)
    testSnp(snpReadHardClip)
  }

  sparkTest("find one snp in reads with 1bp sequence mismatch") {
    val snpRdd = sc.parallelize(Seq(snpReadMCigar,
      snpReadEqCigar,
      snpReadSoftClip,
      snpReadHardClip))
    val variants = DiscoverVariants.variantsInRdd(snpRdd).collect
    assert(variants.size === 1)
    validateSnp(variants.head)
  }

  def validateInsertion(ins: Variant) {
    assert(ins.getContigName() === "2")
    assert(ins.getStart() === 13L)
    assert(ins.getEnd() === 14L)
    assert(ins.getReferenceAllele() === "C")
    assert(ins.getAlternateAllele() === "CTT")
  }

  test("find insertion in read") {
    val variants = DiscoverVariants.variantsInRead(insertRead, 0)
    assert(variants.size === 1)
    validateInsertion(variants.head.toVariant)
  }

  sparkTest("find insertion in reads") {
    val insRdd = sc.parallelize(Seq(insertRead, insertRead, insertRead))
    val variants = DiscoverVariants.variantsInRdd(insRdd).collect
    assert(variants.size === 1)
    validateInsertion(variants.head)
  }

  def validateDeletion(del: Variant) {
    assert(del.getContigName() === "3")
    assert(del.getStart() === 13L)
    assert(del.getEnd() === 16L)
    assert(del.getReferenceAllele() === "CTT")
    assert(del.getAlternateAllele() === "C")
  }

  test("find deletion in read") {
    val variants = DiscoverVariants.variantsInRead(deleteRead, 0)
    assert(variants.size === 1)
    validateDeletion(variants.head.toVariant)
  }

  sparkTest("find deletion in reads") {
    val delRdd = sc.parallelize(Seq(deleteRead, deleteRead, deleteRead))
    val variants = DiscoverVariants.variantsInRdd(delRdd).collect
    assert(variants.size === 1)
    validateDeletion(variants.head)
  }

  sparkTest("find variants in alignment record rdd") {
    val rdd = sc.parallelize(Seq(
      unalignedRead,
      perfectReadMCigar, perfectReadEqCigar,
      snpReadMCigar, snpReadEqCigar,
      insertRead,
      deleteRead))
    val readRdd = AlignmentRecordRDD(rdd,
      SequenceDictionary(
        SequenceRecord("1", 50L),
        SequenceRecord("2", 40L),
        SequenceRecord("3", 30L)),
      RecordGroupDictionary.empty)

    val variantRdd = DiscoverVariants(readRdd)

    assert(variantRdd.rdd.count === 3)
    assert(variantRdd.sequences.records.size === 3)
  }

  test("break TT->CA mnp into two snps") {
    val variants = DiscoverVariants.variantsInRead(mnpRead, 0)
    assert(variants.size === 2)
    assert(variants.forall(_.contigName == "3"))
    assert(variants.forall(_.referenceAllele == "T"))
    val optC = variants.find(_.alternateAllele == "C")
    assert(optC.isDefined)
    optC.foreach(c => {
      assert(c.start === 13)
      assert(c.end === 14)
    })
    val optA = variants.find(_.alternateAllele == "A")
    assert(optA.isDefined)
    optA.foreach(a => {
      assert(a.start === 14)
      assert(a.end === 15)
    })
  }
}
