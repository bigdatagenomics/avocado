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

import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.AlignmentRecord
import scala.collection.mutable.ListBuffer

class HardLimiterSuite extends AvocadoFunSuite {

  val reads = (0 to 5).map(i => {
    AlignmentRecord.newBuilder()
      .setContigName("ctg")
      .setStart(i.toLong)
      .setEnd(i.toLong + 3L)
      .build()
  }).map(r => (r, 1))
    .toSeq

  def addedRead(idx: Int,
                msg: String = "added"): (AlignmentRecord, Int) = {
    val (read, v) = reads(idx)
    (AlignmentRecord.newBuilder(read)
      .setReadName(msg)
      .build(), v)
  }

  test("add a read to an empty buffer") {
    val (emitted, buffer, size) = HardLimiter.processRead(
      reads(0),
      ListBuffer.empty[(AlignmentRecord, Int)],
      0,
      5)

    assert(emitted.isEmpty)
    assert(buffer.size === 1)
    assert(buffer.head === reads(0))
    assert(size === 1)
  }

  test("add a read to a non-empty buffer, without moving forward") {
    val (emitted, buffer, size) = HardLimiter.processRead(
      addedRead(0),
      ListBuffer(reads(0)),
      1,
      5)

    assert(emitted.isEmpty)
    assert(buffer.size === 2)
    assert(buffer(0) === reads(0))
    assert(buffer(1) === addedRead(0))
    assert(size === 2)
  }

  test("add a read to a non-empty buffer, and move forward") {
    val (emitted, buffer, size) = HardLimiter.processRead(
      reads(3),
      ListBuffer(reads(0)),
      1,
      5)

    assert(emitted.size === 1)
    assert(emitted.head === reads(0))
    assert(buffer.size === 1)
    assert(buffer.head === reads(3))
    assert(size === 1)
  }

  test("trying to add a read to a full buffer—without moving forward—drops the read") {
    val (emitted, buffer, size) = HardLimiter.processRead(
      addedRead(2),
      ListBuffer(reads(0), reads(1), reads(2)),
      3,
      3)

    assert(emitted.isEmpty)
    assert(buffer.size === 3)
    assert(buffer(0) === reads(0))
    assert(buffer(1) === reads(1))
    assert(buffer(2) === reads(2))
    assert(size === 3)
  }

  test("add a read to a full buffer, while moving forward and keeping buffer full") {
    val (emitted, buffer, size) = HardLimiter.processRead(
      addedRead(3),
      ListBuffer(reads(0), reads(1), reads(2)),
      3,
      3)

    assert(emitted.size === 1)
    assert(emitted.head === reads(0))
    assert(buffer.size === 3)
    assert(buffer(0) === reads(1))
    assert(buffer(1) === reads(2))
    assert(buffer(2) === addedRead(3))
    assert(size === 3)
  }

  test("add a read to a full buffer, while moving forward and emptying buffer") {
    val (emitted, buffer, size) = HardLimiter.processRead(
      addedRead(5),
      ListBuffer(reads(0), reads(1), reads(2)),
      3,
      3)

    assert(emitted.size === 3)
    assert(emitted(0) === reads(0))
    assert(emitted(1) === reads(1))
    assert(emitted(2) === reads(2))
    assert(buffer.size === 1)
    assert(buffer.head === addedRead(5))
    assert(size === 1)
  }

  test("adding an out of order read should fire an assert") {
    intercept[AssertionError] {
      HardLimiter.processRead(
        reads(0),
        ListBuffer(reads(1)),
        1,
        5)
    }
  }

  test("adding a read that is on the wrong contig should fire an assert") {
    intercept[AssertionError] {
      val randomRead = AlignmentRecord.newBuilder()
        .setContigName("random")
        .setStart(100L)
        .setEnd(101L)
        .build()

      HardLimiter.processRead(
        (randomRead, 1),
        ListBuffer(reads(0)),
        1,
        5)
    }
  }

  test("apply hard limiting to an iterator that is wholly under the coverage limit") {
    val limitedIterator = HardLimiter.limitPartition(reads.toIterator,
      3)
    val iterSize = limitedIterator.size
    assert(iterSize === 6)
  }

  test("apply hard limiting to an iterator that is partially under the coverage limit") {
    val limitedReads = HardLimiter.limitPartition(reads.toIterator,
      2).map(_._1)
      .toSeq

    // we drop every third read
    assert(limitedReads.size === 4)
    val readStarts = limitedReads.map(_.getStart.toInt)
      .toSet
    assert(readStarts(0))
    assert(readStarts(1))
    assert(readStarts(3))
    assert(readStarts(4))
  }

  test("apply hard limiting to an iterator that is wholly over the coverage limit") {
    val readIterator = ((0 to 10).map(i => addedRead(0, msg = i.toString)) ++
      (0 to 6).map(i => addedRead(3, msg = i.toString)).reverse)
      .toIterator
    val reads = HardLimiter.limitPartition(readIterator,
      4).map(_._1).toSeq
    assert(reads.size === 8)
    val (pos0, pos3) = reads.partition(_.getStart < 3L)
    assert(pos0.size === 4)
    assert(pos0(0).getReadName === "0")
    assert(pos0(1).getReadName === "1")
    assert(pos0(2).getReadName === "2")
    assert(pos0(3).getReadName === "3")
    assert(pos3.size === 4)
    assert(pos3(0).getReadName === "6")
    assert(pos3(1).getReadName === "5")
    assert(pos3(2).getReadName === "4")
    assert(pos3(3).getReadName === "3")
  }

  sparkTest("apply hard limiting on a file that is wholly under the coverage limit") {
    val readPath = resourceUrl("NA12878.chr1.832736.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .sortReadsByReferencePosition()

    val limitedReads = reads.transform(rdd => {
      HardLimiter(rdd.map(r => (r, 1)), 100).map(_._1)
    })

    assert(reads.rdd.count === limitedReads.rdd.count)
  }

  sparkTest("apply hard limiting on a file with sections over the coverage limit") {
    val readPath = resourceUrl("NA12878.chr1.832736.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .sortReadsByReferencePosition()

    val limitedReads = reads.transform(rdd => {
      HardLimiter(rdd.map(r => (r, 1)), 50).map(_._1)
    })

    // we should lose four reads
    assert(reads.rdd.count === limitedReads.rdd.count + 4L)

    // no base should be covered by more reads than the threshold
    val coverage = limitedReads.toCoverage()
      .rdd
      .map(_.count)
      .collect
    coverage.forall(_ < 50.0)
  }
}
