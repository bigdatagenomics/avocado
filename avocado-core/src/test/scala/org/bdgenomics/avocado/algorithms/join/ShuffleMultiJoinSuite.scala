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
package org.bdgenomics.avocado.algorithms.join

import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.avocado.AvocadoFunSuite

class ShuffleMultiJoinSuite extends AvocadoFunSuite {

  test("process a partition that starts with a T type") {
    val iterT = Iterator((ReferenceRegion("chr", 0L, 10L), 0),
      (ReferenceRegion("chr", 10L, 20L), 3),
      (ReferenceRegion("chr", 20L, 30L), 6))
    val iterU = Iterator((ReferenceRegion("chr", 0L, 5L), 1L),
      (ReferenceRegion("chr", 5L, 10L), 2L),
      (ReferenceRegion("chr", 10L, 15L), 4L),
      (ReferenceRegion("chr", 15L, 20L), 5L),
      (ReferenceRegion("chr", 20L, 25L), 7L),
      (ReferenceRegion("chr", 25L, 30L), 8L))
    val j = ShuffleMultiJoin.processPartition(iterT, iterU).toSeq

    assert(j.length === 3)
    val kSet = j.map(_._1).toSet
    assert(kSet.size === 3)
    assert(kSet(0))
    assert(kSet(3))
    assert(kSet(6))
    j.foreach(kv => {
      assert(kv._2.size === 2)
      val vSet = kv._2.map(_.toInt).toSet
      assert(vSet(kv._1 + 1))
      assert(vSet(kv._1 + 2))
    })
  }

  test("process a partition that starts with a U type") {
    val iterT = Iterator((ReferenceRegion("chr", 10L, 20L), 3),
      (ReferenceRegion("chr", 20L, 30L), 6))
    val iterU = Iterator((ReferenceRegion("chr", 0L, 5L), 1L),
      (ReferenceRegion("chr", 5L, 10L), 2L),
      (ReferenceRegion("chr", 10L, 15L), 4L),
      (ReferenceRegion("chr", 15L, 20L), 5L),
      (ReferenceRegion("chr", 20L, 25L), 7L),
      (ReferenceRegion("chr", 25L, 30L), 8L))
    val j = ShuffleMultiJoin.processPartition(iterT, iterU).toSeq

    assert(j.length === 2)
    val kSet = j.map(_._1).toSet
    assert(kSet.size === 2)
    assert(kSet(3))
    assert(kSet(6))
    j.foreach(kv => {
      assert(kv._2.size === 2)
      val idx = kv._1
      val vSet = kv._2.map(_.toInt).toSet
      assert(vSet(idx + 1))
      assert(vSet(idx + 2))
    })
  }
}
