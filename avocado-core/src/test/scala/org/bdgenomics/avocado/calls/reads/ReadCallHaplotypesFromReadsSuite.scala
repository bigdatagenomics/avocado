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
package org.bdgenomics.avocado.calls.reads

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichADAMRecord
import org.bdgenomics.adam.util.SparkFunSuite
import org.bdgenomics.formats.avro.{ ADAMGenotypeAllele, ADAMRecord, ADAMContig }
import org.bdgenomics.avocado.algorithms.hmm._
import org.bdgenomics.avocado.partitioners.PartitionSet
import parquet.filter.UnboundRecordFilter
import scala.collection.JavaConversions._
import scala.collection.immutable.SortedMap

class ReadCallHaplotypesFromReadsSuite extends ReadCallHaplotypesSuite {

  val rc_short = new ReadCallHaplotypesFromReads(emptyPartition, 4)
  val rc_long = new ReadCallHaplotypesFromReads(emptyPartition, 20)

  test("Test log sum for similar values") {
    val sum = HaplotypePair.exactLogSumExp10(1.0, 1.0)
    assert(1.3 * 0.99 < sum && 1.3 * 1.01 > sum)
  }

  test("Test log sum for dissimilar values") {
    val sum = HaplotypePair.exactLogSumExp10(1.0, -3.0)
    assert(1.00004342 * 0.99 < sum && 1.00004342 * 1.01 > sum)
  }

  test("insert reads into reference") {
    val reference = "TACCAATGTAA"

    val mismatchReads = Seq(
      make_read("TACCCAT", 0L, "7M", "4A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0),
      make_read("ACCCATG", 1L, "7M", "3A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1),
      make_read("CCCATGT", 2L, "7M", "2A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2))
      .map(r => RichADAMRecord(r))

    val mismatchSeq = mismatchReads.map(rc_short.insertReadIntoReference(_, reference, 0, reference.length))
      .distinct

    assert(mismatchSeq.length === 1)
    mismatchSeq.foreach(s => assert(s === "TACCCATGTAA"))

    val matchReads = Seq(
      make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5),
      make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6),
      make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7))
      .map(r => RichADAMRecord(r))

    val matchSeq = matchReads.map(rc_short.insertReadIntoReference(_, reference, 0, reference.length))
      .distinct

    assert(matchSeq.length === 1)
    matchSeq.foreach(s => assert(s === "TACCAATGTAA"))

    val insertReads = Seq(
      make_read("TACCAAA", 0L, "4M1I2M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0),
      make_read("ACCAAAT", 1L, "3M1I3M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1),
      make_read("CCAAATG", 2L, "2M1I4M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2))
      .map(r => RichADAMRecord(r))

    val insertSeq = insertReads.map(rc_short.insertReadIntoReference(_, reference, 0, reference.length))
      .distinct

    assert(insertSeq.length === 1)
    insertSeq.foreach(s => assert(s === "TACCAAATGTAA"))

    val deleteReads = Seq(
      make_read("TACCATG", 0L, "4M1D3M", "4^A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0),
      make_read("ACCATGT", 1L, "3M1D4M", "3^A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1),
      make_read("CCATGTA", 2L, "2M1D5M", "2^A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2))
      .map(r => RichADAMRecord(r))

    val deleteSeq = deleteReads.map(r => rc_short.insertReadIntoReference(r, reference, 0, reference.length))
      .distinct

    assert(deleteSeq.length === 1)
    deleteSeq.foreach(s => assert(s === "TACCATGTAA"))
  }
}
