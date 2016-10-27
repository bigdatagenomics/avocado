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

import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.{ AlignmentRecord, Variant }

class TreeRegionJoinSuite extends AvocadoFunSuite {

  test("build a forest with a single item and retrieve data") {
    val forest = Forest(Array((ReferenceRegion("chr1", 10L, 15L), 1)))

    assert(forest.length === 1)
    assert(forest.midpoint === 1)

    // retrieve a value wholly inside the first key
    val wholly = forest.get(ReferenceRegion("chr1", 11L, 12L))
    assert(wholly.size === 1)
    assert(wholly.head === 1)

    // retrieve a value that envelops the first key
    val envelops = forest.get(ReferenceRegion("chr1", 5L, 20L))
    assert(envelops.size === 1)
    assert(envelops.head === 1)

    // retrieve a value overlapping the start of the first key
    val start = forest.get(ReferenceRegion("chr1", 7L, 11L))
    assert(start.size === 1)
    assert(start.head === 1)

    // retrieve a value overlapping the end of the first key
    val end = forest.get(ReferenceRegion("chr1", 14L, 16L))
    assert(end.size === 1)
    assert(end.head === 1)

    // retrieve a value before the first key
    val before = forest.get(ReferenceRegion("chr1", 2L, 5L))
    assert(before.isEmpty)

    // retrieve a value after the first key
    val after = forest.get(ReferenceRegion("chr1", 22L, 75L))
    assert(after.isEmpty)

    // retrieve a value on a different contig
    val otherContig = forest.get(ReferenceRegion("chr5", 10L, 14L))
    assert(otherContig.isEmpty)
  }

  test("build a forest with data from a single contig and retrieve data") {
    val forest = Forest(Array((ReferenceRegion("chr1", 10L, 15L), 0),
      (ReferenceRegion("chr1", 12L, 20L), 1),
      (ReferenceRegion("chr1", 24L, 30L), 2),
      (ReferenceRegion("chr1", 55L, 65L), 3)))

    assert(forest.length === 4)
    assert(forest.midpoint === 2)

    // retrieve a value solely inside the first key
    val solely = forest.get(ReferenceRegion("chr1", 11L, 12L))
    assert(solely.size === 1)
    assert(solely.head === 0)

    // retrieve a value wholly inside the first and second keys
    val wholly = forest.get(ReferenceRegion("chr1", 12L, 13L)).toSet
    assert(wholly.size === 2)
    assert(wholly(0))
    assert(wholly(1))

    // retrieve a value that envelops the third and fourth keys
    val envelops = forest.get(ReferenceRegion("chr1", 20L, 100L)).toSet
    assert(envelops.size === 2)
    assert(envelops(2))
    assert(envelops(3))

    // retrieve a value overlapping the start of the fourth key
    val start = forest.get(ReferenceRegion("chr1", 50L, 60L))
    assert(start.size === 1)
    assert(start.head === 3)

    // retrieve a value overlapping the end of the third key
    val end = forest.get(ReferenceRegion("chr1", 26L, 36L))
    assert(end.size === 1)
    assert(end.head === 2)

    // retrieve a value before the first key
    val before = forest.get(ReferenceRegion("chr1", 2L, 5L))
    assert(before.isEmpty)

    // retrieve a value between the second and third keys
    val between = forest.get(ReferenceRegion("chr1", 21L, 22L))
    assert(between.isEmpty)

    // retrieve a value after the last key
    val after = forest.get(ReferenceRegion("chr1", 70L, 75L))
    assert(after.isEmpty)

    // retrieve a value on a different contig
    val otherContig = forest.get(ReferenceRegion("chr5", 10L, 14L))
    assert(otherContig.isEmpty)
  }

  test("build a forest with data from multiple contigs and retrieve data") {
    val forest = Forest(Array((ReferenceRegion("chr1", 10L, 15L), 0),
      (ReferenceRegion("chr1", 12L, 20L), 1),
      (ReferenceRegion("chr2", 24L, 30L), 2),
      (ReferenceRegion("chr2", 55L, 65L), 3),
      (ReferenceRegion("chr2", 75L, 80L), 4),
      (ReferenceRegion("chr3", 10L, 15L), 5)))

    assert(forest.length === 6)
    assert(forest.midpoint === 4)

    // get a value that overlaps just the second key
    val second = forest.get(ReferenceRegion("chr1", 15L, 25L))
    assert(second.size === 1)
    assert(second.head === 1)

    // get all values on the second contig
    val secondContig = forest.get(ReferenceRegion("chr2", 25L, 80L)).toSet
    assert(secondContig.size === 3)
    assert(secondContig(2))
    assert(secondContig(3))
    assert(secondContig(4))

    // get the value on the last contig
    val last = forest.get(ReferenceRegion("chr3", 5L, 12L))
    assert(last.size === 1)
    assert(last.head === 5)
  }

  sparkTest("build a forest out of data on a single contig and retrieve data") {
    val rdd = sc.parallelize(Seq((ReferenceRegion("chr1", 10L, 15L), 1),
      (ReferenceRegion("chr1", 9L, 12L), 0),
      (ReferenceRegion("chr1", 100L, 150L), 4),
      (ReferenceRegion("chr1", 80L, 95L), 2),
      (ReferenceRegion("chr1", 80L, 110L), 3)))

    val forest = Forest(rdd)

    assert(forest.length === 5)
    assert(forest.midpoint === 4)
    (0 until forest.length).foreach(idx => {
      assert(forest.array(idx)._2 === idx)
    })

    // retrieve a value overlapping the first two keys
    val firstTwo = forest.get(ReferenceRegion("chr1", 10L, 12L)).toSet
    assert(firstTwo.size === 2)
    assert(firstTwo(0))
    assert(firstTwo(1))

    // retrieve a value overlapping the last three keys
    val lastThree = forest.get(ReferenceRegion("chr1", 90L, 105L)).toSet
    assert(lastThree.size === 3)
    assert(lastThree(2))
    assert(lastThree(3))
    assert(lastThree(4))

    // retrieve a value overlapping just the last key
    val last = forest.get(ReferenceRegion("chr1", 130L, 140L))
    assert(last.size === 1)
    assert(last.head === 4)

    // retrieve a value before the first key
    val before = forest.get(ReferenceRegion("chr1", 2L, 5L))
    assert(before.isEmpty)

    // retrieve a value between the second and third keys
    val between = forest.get(ReferenceRegion("chr1", 21L, 22L))
    assert(between.isEmpty)

    // retrieve a value after the last key
    val after = forest.get(ReferenceRegion("chr1", 500L, 675L))
    assert(after.isEmpty)

    // retrieve a value on a different contig
    val otherContig = forest.get(ReferenceRegion("chr5", 10L, 14L))
    assert(otherContig.isEmpty)
  }

  sparkTest("run a join between data on a single contig") {

    val rightRdd = sc.parallelize(Seq(
      (ReferenceRegion("chr1", 10L, 20L), 0),
      (ReferenceRegion("chr1", 15L, 25L), 1),
      (ReferenceRegion("chr1", 30L, 50L), 2),
      (ReferenceRegion("chr1", 60L, 70L), 3),
      (ReferenceRegion("chr1", 90L, 100L), 4)))
      .map(kv => {
        val (k, v) = kv
        // i have made many poor life decisions
        (k, Variant.newBuilder
          .setStart(v.toLong)
          .build)
      })

    val leftRdd = sc.parallelize(Seq(
      (ReferenceRegion("chr1", 12L, 22L), 0),
      (ReferenceRegion("chr1", 20L, 35L), 1),
      (ReferenceRegion("chr1", 40L, 55L), 2),
      (ReferenceRegion("chr1", 75L, 85L), 3),
      (ReferenceRegion("chr1", 95L, 105L), 4)))
      .map(kv => {
        val (k, v) = kv
        // and this is but another one of them
        (k, AlignmentRecord.newBuilder
          .setStart(v.toLong)
          .build)
      })

    val joinData = TreeRegionJoin.joinAndGroupByRight(rightRdd, leftRdd)
      .map(kv => {
        val (k, v) = kv
        (k.map(_.getStart.toInt), v.getStart.toInt)
      }).collect

    assert(joinData.size === 4)

    val joinMap = joinData.map(_.swap)
      .toMap
      .mapValues(_.toSet)

    assert(joinMap.size === 4)
    assert(joinMap(0).size === 2)
    assert(joinMap(0)(0))
    assert(joinMap(0)(1))
    assert(joinMap(1).size === 2)
    assert(joinMap(1)(1))
    assert(joinMap(1)(2))
    assert(joinMap(2).size === 1)
    assert(joinMap(2)(2))
    assert(joinMap(4).size === 1)
    assert(joinMap(4)(4))
  }
}
