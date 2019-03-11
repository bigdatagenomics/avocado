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
package org.bdgenomics.avocado.models

import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.feature.FeatureDataset
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.Feature

class CopyNumberMapSuite extends AvocadoFunSuite {

  test("create an empty map") {
    val emptyMap = CopyNumberMap.empty(2)

    assert(emptyMap.basePloidy === 2)
    assert(emptyMap.minPloidy === 2)
    assert(emptyMap.maxPloidy === 2)
    assert(emptyMap.variantsByReference.isEmpty)
  }

  sparkTest("create a map with only diploid features") {
    val cnvs = Seq(Feature.newBuilder
      .setStart(100L)
      .setEnd(201L)
      .setReferenceName("chr1")
      .setFeatureType("DIP")
      .build)

    val emptyMap = CopyNumberMap(2,
      FeatureDataset(sc.parallelize(cnvs)))

    assert(emptyMap.basePloidy === 2)
    assert(emptyMap.minPloidy === 2)
    assert(emptyMap.maxPloidy === 2)
    assert(emptyMap.variantsByReference.isEmpty)
    assert(emptyMap.overlappingVariants(ReferenceRegion("chr1", 100L, 201L))
      .isEmpty)
  }

  sparkTest("create a map with a mix of features") {
    val cnvs = Seq(Feature.newBuilder
      .setStart(100L)
      .setEnd(201L)
      .setReferenceName("chr1")
      .setFeatureType("DIP")
      .build,
      Feature.newBuilder
        .setStart(1000L)
        .setEnd(2000L)
        .setReferenceName("chr1")
        .setFeatureType("DUP")
        .build,
      Feature.newBuilder
        .setStart(2000L)
        .setEnd(3000L)
        .setReferenceName("chr1")
        .setFeatureType("DEL")
        .build,
      Feature.newBuilder
        .setStart(2000L)
        .setEnd(3000L)
        .setReferenceName("chr2")
        .setFeatureType("DEL")
        .build)

    val cnvMap = CopyNumberMap(2,
      FeatureDataset(sc.parallelize(cnvs)))

    assert(cnvMap.basePloidy === 2)
    assert(cnvMap.minPloidy === 1)
    assert(cnvMap.maxPloidy === 3)
    assert(cnvMap.variantsByReference.size === 2)
    val chr1Cnvs = cnvMap.variantsByReference("chr1")
    assert(chr1Cnvs.size === 2)
    assert(chr1Cnvs(0)._1 === ReferenceRegion("chr1", 1000L, 2000L))
    assert(chr1Cnvs(0)._2 === 3)
    assert(chr1Cnvs(1)._1 === ReferenceRegion("chr1", 2000L, 3000L))
    assert(chr1Cnvs(1)._2 === 1)
    val chr2Cnvs = cnvMap.variantsByReference("chr2")
    assert(chr2Cnvs.size === 1)
    assert(chr2Cnvs(0)._1 === ReferenceRegion("chr2", 2000L, 3000L))
    assert(chr2Cnvs(0)._2 === 1)

    assert(cnvMap.overlappingVariants(ReferenceRegion("chr1", 100L, 150L))
      .isEmpty)
    assert(cnvMap.overlappingVariants(ReferenceRegion("chr1", 100L, 2000L))
      .size === 1)
    assert(cnvMap.overlappingVariants(ReferenceRegion("chr1", 100L, 2000L))
      .head._2 === 3)
    val allChr1Cnvs = cnvMap.overlappingVariants(ReferenceRegion("chr1",
      100L,
      2001L)).map(_._2).toSet
    assert(allChr1Cnvs.size === 2)
    assert(allChr1Cnvs(1))
    assert(allChr1Cnvs(3))
    assert(cnvMap.overlappingVariants(ReferenceRegion("chr1", 2000L, 2001L))
      .size === 1)
    assert(cnvMap.overlappingVariants(ReferenceRegion("chr1", 2000L, 2001L))
      .head._2 === 1)
    assert(cnvMap.overlappingVariants(ReferenceRegion("chr2", 100L, 2001L))
      .size === 1)
    assert(cnvMap.overlappingVariants(ReferenceRegion("chr2", 100L, 2001L))
      .head._2 === 1)
  }
}
