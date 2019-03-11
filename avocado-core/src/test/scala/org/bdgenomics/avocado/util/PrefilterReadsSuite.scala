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
package org.bdgenomics.avocado.util

import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.adam.models.{
  ReadGroupDictionary,
  SequenceDictionary,
  SequenceRecord
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.AlignmentRecordDataset
import org.bdgenomics.formats.avro.AlignmentRecord

case class TestPrefilterReadsArgs(var autosomalOnly: Boolean = false,
                                  var keepMitochondrialChromosome: Boolean = false,
                                  var keepDuplicates: Boolean = true,
                                  var minMappingQuality: Int = -1,
                                  var keepNonPrimary: Boolean = true) extends PrefilterReadsArgs {
}

class PrefilterReadsSuite extends AvocadoFunSuite {

  test("filter on read uniqueness") {
    val uniqueRead = AlignmentRecord.newBuilder
      .setDuplicateRead(false)
      .build
    assert(PrefilterReads.filterUnique(uniqueRead))

    val duplicateRead = AlignmentRecord.newBuilder
      .setDuplicateRead(true)
      .build
    assert(!PrefilterReads.filterUnique(duplicateRead))
  }

  test("filter unmapped reads") {
    val mappedRead = AlignmentRecord.newBuilder
      .setReadMapped(true)
      .setPrimaryAlignment(true)
      .build
    assert(PrefilterReads.filterMapped(mappedRead, false))

    val nonPrimaryMappedRead = AlignmentRecord.newBuilder
      .setReadMapped(true)
      .setPrimaryAlignment(false)
      .build
    assert(!PrefilterReads.filterMapped(nonPrimaryMappedRead, false))
    assert(PrefilterReads.filterMapped(nonPrimaryMappedRead, true))

    val unmappedRead = AlignmentRecord.newBuilder
      .setReadMapped(false)
      .build
    assert(!PrefilterReads.filterMapped(unmappedRead, true))
  }

  val referenceNames = Seq("chr1",
    "1",
    "chrX",
    "X",
    "chrY",
    "Y",
    "chrM",
    "MT")

  def testChromosomeHelper(testFn: (String => Boolean), passIdx: Int) {
    testChromosomeHelperSet(testFn, Set(passIdx))
  }

  def testChromosomeHelperSet(testFn: (String => Boolean), passIdxSet: Set[Int]) {
    def assertIdx(idx: Int, testStr: String) = {
      if (passIdxSet(idx)) {
        assert(testFn(testStr))
      } else {
        assert(!testFn(testStr))
      }
    }

    referenceNames.zipWithIndex
      .foreach(p => assertIdx(p._2, p._1))
  }

  test("filter autosomal chromosomes with grc names") {
    testChromosomeHelper(PrefilterReads.filterGrcAutosome, 0)
  }

  test("filter sex chromosomes with grc names") {
    testChromosomeHelperSet(PrefilterReads.filterGrcSex, Set(2, 4))
  }

  test("filter mitochondrial chromosome with a grc names") {
    testChromosomeHelper(PrefilterReads.filterGrcMitochondrial, 6)
  }

  test("filter autosomal chromosomes with hg names") {
    testChromosomeHelper(PrefilterReads.filterNonGrcAutosome, 1)
  }

  test("filter sex chromosomes with hg names") {
    testChromosomeHelperSet(PrefilterReads.filterNonGrcSex, Set(3, 5))
  }

  test("filter mitochondrial chromosome with a hg names") {
    testChromosomeHelper(PrefilterReads.filterNonGrcMitochondrial, 7)
  }

  test("filter autosomal chromosomes from generator") {
    testChromosomeHelperSet(PrefilterReads.referenceFilterFn(TestPrefilterReadsArgs(autosomalOnly = true)), Set(0, 1))
  }

  test("filter autosomal + sex chromosomes from generator") {
    testChromosomeHelperSet(PrefilterReads.referenceFilterFn(TestPrefilterReadsArgs()), Set(0, 1,
      2, 3,
      4, 5))
  }

  test("filter all chromosomes from generator") {
    testChromosomeHelperSet(PrefilterReads.referenceFilterFn(TestPrefilterReadsArgs(keepMitochondrialChromosome = true)), Set(0, 1, 2, 3, 4, 5, 6, 7))
  }

  test("update a read whose mate is mapped to a filtered reference") {
    val read = AlignmentRecord.newBuilder()
      .setReadPaired(true)
      .setMateMapped(true)
      .setMateReferenceName("notARealReference")
      .build

    val filters = PrefilterReads.referenceFilterFn(TestPrefilterReadsArgs())
    val nullified = PrefilterReads.maybeNullifyMate(read, filters)

    assert(!nullified.getMateMapped)
    assert(nullified.getMateReferenceName == null)
  }

  val reads = Seq(AlignmentRecord.newBuilder()
    .setReadMapped(false),
    AlignmentRecord.newBuilder()
      .setReadMapped(true)
      .setDuplicateRead(true),
    AlignmentRecord.newBuilder()
      .setReadMapped(true)
      .setDuplicateRead(false)).flatMap(rb => {
      referenceNames.map(cn => rb.setReferenceName(cn).build)
    })

  def testReadHelperSet(testArgs: PrefilterReadsArgs, passIdxSet: Set[Int]) {
    val testFn = PrefilterReads.readFilterFn(testArgs,
      PrefilterReads.referenceFilterFn(testArgs))

    def assertIdx(idx: Int, testRead: AlignmentRecord) = {
      if (passIdxSet(idx)) {
        assert(testFn(testRead))
      } else {
        assert(!testFn(testRead))
      }
    }

    reads.zipWithIndex
      .foreach(p => assertIdx(p._2, p._1))
  }

  test("filter reads mapped to autosomal chromosomes from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(autosomalOnly = true), Set(8, 9, 16, 17))
  }

  test("filter reads mapped to autosomal + sex chromosomes from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(), Set(8, 9, 10, 11, 12, 13,
      16, 17, 18, 19, 20, 21))
  }

  test("filter reads mapped to all chromosomes from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepMitochondrialChromosome = true),
      Set(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23))
  }

  test("filter reads uniquely mapped to autosomal chromosomes from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false,
      autosomalOnly = true), Set(16, 17))
  }

  test("filter reads uniquely mapped to autosomal + sex chromosomes from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false),
      Set(16, 17, 18, 19, 20, 21))
  }

  test("filter reads uniquely mapped to all chromosomes from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false,
      keepMitochondrialChromosome = true),
      Set(16, 17, 18, 19, 20, 21, 22, 23))
  }

  val sequences = new SequenceDictionary(referenceNames.map(cn => SequenceRecord(cn, 10L))
    .toVector)

  def testRdd(args: PrefilterReadsArgs, numReads: Int, numReferences: Int) {

    val readRdd = AlignmentRecordDataset(sc.parallelize(reads), sequences, ReadGroupDictionary.empty, Seq.empty)
    val filteredRdd = PrefilterReads(readRdd, args)

    assert(filteredRdd.rdd.count === numReads)
    assert(filteredRdd.sequences.records.size === numReferences)
  }

  sparkTest("filter rdd of reads mapped to autosomal chromosomes from generator") {
    testRdd(TestPrefilterReadsArgs(autosomalOnly = true), 4, 2)
  }

  sparkTest("filter rdd of reads mapped to autosomal + sex chromosomes from generator") {
    testRdd(TestPrefilterReadsArgs(), 12, 6)
  }

  sparkTest("filter rdd of reads mapped to all chromosomes from generator") {
    testRdd(TestPrefilterReadsArgs(keepMitochondrialChromosome = true), 16, 8)
  }

  sparkTest("filter rdd of reads uniquely mapped to autosomal chromosomes from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false,
      autosomalOnly = true), 2, 2)
  }

  sparkTest("filter rdd of reads uniquely mapped to autosomal + sex chromosomes from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false),
      6, 6)
  }

  sparkTest("filter rdd of reads uniquely mapped to all chromosomes from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false,
      keepMitochondrialChromosome = true),
      8, 8)
  }
}
