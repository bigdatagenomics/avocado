/*
 * Copyright (c) 2014. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.avocado.partitioners

import org.bdgenomics.adam.models.ReferenceRegion
import org.scalatest.FunSuite
import scala.collection.immutable.SortedMap

class DefaultPartitionerSuite extends FunSuite {

  val contigs = Map(0 -> 1000L, 1 -> 500L, 2 -> 100L)
  val partitioner = new DefaultPartitioner(contigs, 100)
  val partitions = partitioner.computePartitions()

  test("build a single partiton per contig") {
    val regions = contigs.toList.map(kv => ReferenceRegion(kv._1, 0, kv._2))

    regions.foreach(r => {
      assert(partitions.isInSet(r))
      assert(!partitions.isOutsideOfSet(r))
    })
  }

  test("check that offsets are correct") {
    val ctg0 = ReferenceRegion(0, 0, 99L)
    assert(partitions.isInSet(ctg0))
    assert(!partitions.isOutsideOfSet(ctg0))
    assert(partitions.getPartition(ctg0).length === 1)
    assert(partitions.getPartition(ctg0).head === 0)

    val ctg1 = ReferenceRegion(1, 0, 99L)
    assert(partitions.isInSet(ctg1))
    assert(!partitions.isOutsideOfSet(ctg1))
    assert(partitions.getPartition(ctg1).length === 1)
    assert(partitions.getPartition(ctg1).head === 10)

    val ctg2 = ReferenceRegion(2, 0, 99L)
    assert(partitions.isInSet(ctg2))
    assert(!partitions.isOutsideOfSet(ctg2))
    assert(partitions.getPartition(ctg2).length === 1)
    assert(partitions.getPartition(ctg2).head === 15)
  }

  test("default partition set properly handles regions that span buckets") {
    val ctg0 = ReferenceRegion(0, 0, 299L)
    assert(partitions.isInSet(ctg0))
    assert(!partitions.isOutsideOfSet(ctg0))
    assert(partitions.getPartition(ctg0).length === 3)
    assert(partitions.getPartition(ctg0).contains(0))
    assert(partitions.getPartition(ctg0).contains(1))
    assert(partitions.getPartition(ctg0).contains(2))
  }

  test("throws exception for contigs not in reference") {
    val ctg3 = ReferenceRegion(3, 0, 100L)

    intercept[IllegalArgumentException] {
      partitions.getPartition(ctg3)
    }
  }

  test("throws exception for contigs that extend past ends of reference contigs") {
    val ctg0 = ReferenceRegion(0, 900L, 1100L)
    assert(partitions.isInSet(ctg0))

    intercept[IllegalArgumentException] {
      partitions.getPartition(ctg0)
    }
  }
}
