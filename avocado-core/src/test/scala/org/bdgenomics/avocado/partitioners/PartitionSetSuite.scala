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

class PartitionSetSuite extends FunSuite {

  val regionMapping = SortedMap(ReferenceRegion(0, 0L, 100L) -> 0,
    ReferenceRegion(0, 100L, 200L) -> 1,
    ReferenceRegion(0, 300L, 400L) -> 2)
  val partitionSet = new PartitionSet(regionMapping)

  test("region that is fully contained in a single partition") {
    val r = ReferenceRegion(0, 10L, 40L)

    assert(partitionSet.isInSet(r))
    assert(!partitionSet.isOutsideOfSet(r))
    assert(partitionSet.getPartition(r).length === 1)
    assert(partitionSet.getPartition(r).head === 0)
  }

  test("region that is fully contained in two adjacent partitions") {
    val r = ReferenceRegion(0, 90L, 120L)

    assert(partitionSet.isInSet(r))
    assert(!partitionSet.isOutsideOfSet(r))
    assert(partitionSet.getPartition(r).length === 2)
    assert(partitionSet.getPartition(r).contains(0))
    assert(partitionSet.getPartition(r).contains(1))
  }

  test("region that is in two adjacent partitions, but extends outside of them as well") {
    val r = ReferenceRegion(0, 90L, 220L)

    assert(partitionSet.isInSet(r))
    assert(partitionSet.isOutsideOfSet(r))
    assert(partitionSet.getPartition(r).length === 2)
    assert(partitionSet.getPartition(r).contains(0))
    assert(partitionSet.getPartition(r).contains(1))
  }

  test("region that is in two partitions with a gap between them") {
    val r = ReferenceRegion(0, 170L, 320L)

    assert(partitionSet.isInSet(r))
    assert(partitionSet.isOutsideOfSet(r))
    assert(partitionSet.getPartition(r).length === 2)
    assert(partitionSet.getPartition(r).contains(1))
    assert(partitionSet.getPartition(r).contains(2))
  }

  test("region that is fully outside of the set") {
    val r = ReferenceRegion(0, 450L, 550L)

    assert(!partitionSet.isInSet(r))
    assert(partitionSet.isOutsideOfSet(r))
    assert(partitionSet.getPartition(r).length === 0)
  }
}
