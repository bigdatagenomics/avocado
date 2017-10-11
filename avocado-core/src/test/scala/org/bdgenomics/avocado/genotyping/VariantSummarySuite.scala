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

import org.bdgenomics.formats.avro.{
  Genotype,
  Variant,
  VariantAnnotation
}
import org.scalatest.FunSuite
import scala.collection.JavaConversions._

class VariantSummarySuite extends FunSuite {

  test("create from genotype without strand bias components") {
    val gt = Genotype.newBuilder
      .setReadDepth(10)
      .setReferenceReadDepth(5)
      .build

    val vs = VariantSummary(gt)

    assert(vs.readDepth === Some(10))
    assert(vs.referenceReadDepth === Some(5))
    assert(vs.forwardReadDepth.isEmpty)
    assert(vs.reverseReadDepth.isEmpty)
    assert(vs.forwardReferenceReadDepth.isEmpty)
    assert(vs.reverseReferenceReadDepth.isEmpty)
  }

  test("create from genotype with strand bias components") {
    val gt = Genotype.newBuilder
      .setReadDepth(10)
      .setReferenceReadDepth(5)
      .setStrandBiasComponents(Seq(2, 3, 4, 1).map(i => i: java.lang.Integer))
      .build

    val vs = VariantSummary(gt)

    assert(vs.readDepth === Some(10))
    assert(vs.referenceReadDepth === Some(5))
    assert(vs.forwardReadDepth === Some(6))
    assert(vs.reverseReadDepth === Some(4))
    assert(vs.forwardReferenceReadDepth === Some(2))
    assert(vs.reverseReferenceReadDepth === Some(3))
  }

  test("invalid strand bias causes exception") {
    val gt = Genotype.newBuilder
      .setStrandBiasComponents(Seq(0, 1, 2).map(i => i: java.lang.Integer))
      .build

    intercept[IllegalArgumentException] {
      VariantSummary(gt)
    }
  }

  test("merge two fully populated summaries") {
    val vs1 = VariantSummary(Some(10),
      Some(5),
      Some(6),
      Some(4),
      Some(3),
      Some(2))
    val vs2 = VariantSummary(Some(12),
      Some(7),
      Some(6),
      Some(4),
      Some(6),
      Some(2))

    val va = vs1.merge(vs2)
      .toAnnotation(Variant.newBuilder.build())

    assert(va.getReadDepth === 22)
    assert(va.getReferenceReadDepth === 12)
    assert(va.getForwardReadDepth === 12)
    assert(va.getReverseReadDepth === 8)
    assert(va.getReferenceForwardReadDepth === 9)
    assert(va.getReferenceReverseReadDepth === 4)
  }

  test("merge two partially populated summaries") {
    val vs1 = VariantSummary(Some(0),
      Some(1),
      Some(1),
      None,
      Some(3),
      None)
    val vs2 = VariantSummary(Some(0),
      None,
      Some(1),
      None,
      None,
      None)
    val vm = vs1.merge(vs2)

    assert(vm.readDepth === Some(0))
    assert(vm.referenceReadDepth === Some(1))
    assert(vm.forwardReadDepth === Some(2))
    assert(vm.reverseReadDepth.isEmpty)
    assert(vm.forwardReferenceReadDepth === Some(3))
    assert(vm.reverseReferenceReadDepth.isEmpty)
  }

  test("populating an annotation should carry old fields") {
    val va = VariantAnnotation.newBuilder
      .setReadDepth(6)
      .setReferenceReadDepth(2)
      .setForwardReadDepth(3)
      .setReferenceForwardReadDepth(2)
      .setReverseReadDepth(1)
      .setReferenceReverseReadDepth(0)
      .build
    val v = Variant.newBuilder
      .setAnnotation(va)
      .build

    val newVa = VariantSummary(None, None,
      None, None,
      None, None)
      .toAnnotation(v)

    assert(newVa === va)
  }
}
