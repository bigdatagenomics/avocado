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
  RecordGroupDictionary,
  SequenceDictionary,
  SequenceRecord
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.formats.avro.AlignmentRecord

case class TestPrefilterReadsArgs(var isNotGrc: Boolean = false,
                                  var autosomalOnly: Boolean = false,
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

  val contigNames = Seq("chr1",
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

    contigNames.zipWithIndex
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

  test("filter autosomal chromosomes with grc names from generator") {
    testChromosomeHelper(PrefilterReads.contigFilterFn(TestPrefilterReadsArgs(autosomalOnly = true)), 0)
  }

  test("filter autosomal + sex chromosomes with grc names from generator") {
    testChromosomeHelperSet(PrefilterReads.contigFilterFn(TestPrefilterReadsArgs()), Set(0, 2, 4))
  }

  test("filter all chromosomes with grc names from generator") {
    testChromosomeHelperSet(PrefilterReads.contigFilterFn(TestPrefilterReadsArgs(keepMitochondrialChromosome = true)), Set(0, 2, 4, 6))
  }

  test("filter autosomal chromosomes with hg names from generator") {
    testChromosomeHelper(PrefilterReads.contigFilterFn(TestPrefilterReadsArgs(isNotGrc = true,
      autosomalOnly = true)), 1)
  }

  test("filter autosomal + sex chromosomes with hg names from generator") {
    testChromosomeHelperSet(PrefilterReads.contigFilterFn(TestPrefilterReadsArgs(isNotGrc = true)), Set(1, 3, 5))
  }

  test("filter all chromosomes with hg names from generator") {
    testChromosomeHelperSet(PrefilterReads.contigFilterFn(TestPrefilterReadsArgs(isNotGrc = true,
      keepMitochondrialChromosome = true)), Set(1, 3, 5, 7))
  }

  val reads = Seq(AlignmentRecord.newBuilder()
    .setReadMapped(false),
    AlignmentRecord.newBuilder()
      .setReadMapped(true)
      .setDuplicateRead(true),
    AlignmentRecord.newBuilder()
      .setReadMapped(true)
      .setDuplicateRead(false)).flatMap(rb => {
      contigNames.map(cn => rb.setContigName(cn).build)
    })

  def testReadHelperSet(testArgs: PrefilterReadsArgs, passIdxSet: Set[Int]) {
    val testFn = PrefilterReads.readFilterFn(testArgs,
      PrefilterReads.contigFilterFn(testArgs))

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

  test("filter reads mapped to autosomal chromosomes with grc names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(autosomalOnly = true), Set(8, 16))
  }

  test("filter reads mapped to autosomal + sex chromosomes with grc names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(), Set(8, 10, 12, 16, 18, 20))
  }

  test("filter reads mapped to all chromosomes with grc names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepMitochondrialChromosome = true),
      Set(8, 10, 12, 14, 16, 18, 20, 22))
  }

  test("filter reads mapped to autosomal chromosomes with hg names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(isNotGrc = true,
      autosomalOnly = true), Set(9, 17))
  }

  test("filter reads mapped to autosomal + sex chromosomes with hg names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(isNotGrc = true), Set(9, 11, 13, 17, 19, 21))
  }

  test("filter reads mapped to all chromosomes with hg names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(isNotGrc = true,
      keepMitochondrialChromosome = true),
      Set(9, 11, 13, 15, 17, 19, 21, 23))
  }

  test("filter reads uniquely mapped to autosomal chromosomes with grc names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false,
      autosomalOnly = true), Set(16))
  }

  test("filter reads uniquely mapped to autosomal + sex chromosomes with grc names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false),
      Set(16, 18, 20))
  }

  test("filter reads uniquely mapped to all chromosomes with grc names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false,
      keepMitochondrialChromosome = true),
      Set(16, 18, 20, 22))
  }

  test("filter reads uniquely mapped to autosomal chromosomes with hg names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false,
      isNotGrc = true,
      autosomalOnly = true), Set(17))
  }

  test("filter reads uniquely mapped to autosomal + sex chromosomes with hg names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false,
      isNotGrc = true), Set(17, 19, 21))
  }

  test("filter reads uniquely mapped to all chromosomes with hg names from generator") {
    testReadHelperSet(TestPrefilterReadsArgs(keepDuplicates = false,
      isNotGrc = true,
      keepMitochondrialChromosome = true),
      Set(17, 19, 21, 23))
  }

  val sequences = new SequenceDictionary(contigNames.map(cn => SequenceRecord(cn, 10L))
    .toVector)

  def testRdd(args: PrefilterReadsArgs, numReads: Int, numContigs: Int) {

    val readRdd = AlignmentRecordRDD(sc.parallelize(reads), sequences, RecordGroupDictionary.empty)
    val filteredRdd = PrefilterReads(readRdd, args)

    assert(filteredRdd.rdd.count === numReads)
    assert(filteredRdd.sequences.records.size === numContigs)
  }

  sparkTest("filter rdd of reads mapped to autosomal chromosomes with grc names from generator") {
    testRdd(TestPrefilterReadsArgs(autosomalOnly = true), 2, 1)
  }

  sparkTest("filter rdd of reads mapped to autosomal + sex chromosomes with grc names from generator") {
    testRdd(TestPrefilterReadsArgs(), 6, 3)
  }

  sparkTest("filter rdd of reads mapped to all chromosomes with grc names from generator") {
    testRdd(TestPrefilterReadsArgs(keepMitochondrialChromosome = true), 8, 4)
  }

  sparkTest("filter rdd of reads mapped to autosomal chromosomes with hg names from generator") {
    testRdd(TestPrefilterReadsArgs(isNotGrc = true,
      autosomalOnly = true), 2, 1)
  }

  sparkTest("filter rdd of reads mapped to autosomal + sex chromosomes with hg names from generator") {
    testRdd(TestPrefilterReadsArgs(isNotGrc = true), 6, 3)
  }

  sparkTest("filter rdd of reads mapped to all chromosomes with hg names from generator") {
    testRdd(TestPrefilterReadsArgs(isNotGrc = true,
      keepMitochondrialChromosome = true),
      8, 4)
  }

  sparkTest("filter rdd of reads uniquely mapped to autosomal chromosomes with grc names from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false,
      autosomalOnly = true), 1, 1)
  }

  sparkTest("filter rdd of reads uniquely mapped to autosomal + sex chromosomes with grc names from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false),
      3, 3)
  }

  sparkTest("filter rdd of reads uniquely mapped to all chromosomes with grc names from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false,
      keepMitochondrialChromosome = true),
      4, 4)
  }

  sparkTest("filter rdd of reads uniquely mapped to autosomal chromosomes with hg names from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false,
      isNotGrc = true,
      autosomalOnly = true), 1, 1)
  }

  sparkTest("filter rdd of reads uniquely mapped to autosomal + sex chromosomes with hg names from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false,
      isNotGrc = true), 3, 3)
  }

  sparkTest("filter rdd of reads uniquely mapped to all chromosomes with hg names from generator") {
    testRdd(TestPrefilterReadsArgs(keepDuplicates = false,
      isNotGrc = true,
      keepMitochondrialChromosome = true),
      4, 4)
  }
}
