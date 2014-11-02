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
package org.bdgenomics.avocado.postprocessing

import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.formats.avro.{ Genotype, GenotypeAllele, Variant }
import org.scalatest.FunSuite

class ReferenceCallFilterSuite extends FunSuite {

  val rcf = new ReferenceCallFilter()

  val ev = Variant.newBuilder().build()

  test("filter out hom ref call and no call") {
    val calls = rcf.filterGenotypes(Seq(Genotype.newBuilder()
      .setVariant(ev)
      .setAlleles(Array(GenotypeAllele.Ref, GenotypeAllele.Ref, GenotypeAllele.Ref).toList)
      .build(),
      Genotype.newBuilder()
        .setVariant(ev)
        .setAlleles(Array(GenotypeAllele.NoCall).toList)
        .build()))

    assert(calls.size === 0)
  }

  test("don't filter out hom alt, biallelic, and het calls") {
    assert(rcf.filterGenotypes(Seq(Genotype.newBuilder()
      .setVariant(ev)
      .setAlleles(Array(GenotypeAllele.Ref, GenotypeAllele.Alt).toList)
      .build())).size === 1)
    assert(rcf.filterGenotypes(Seq(Genotype.newBuilder()
      .setVariant(ev)
      .setAlleles(Array(GenotypeAllele.Alt, GenotypeAllele.OtherAlt).toList)
      .build())).size === 1)
    assert(rcf.filterGenotypes(Seq(Genotype.newBuilder()
      .setVariant(ev)
      .setAlleles(Array(GenotypeAllele.Alt).toList)
      .build())).size === 1)
  }

  test("don't filter ref calls at a site with non-ref calls") {
    val calls = rcf.filterGenotypes(Seq(Genotype.newBuilder()
      .setVariant(ev)
      .setAlleles(Array(GenotypeAllele.Ref, GenotypeAllele.Ref).toList)
      .build(),
      Genotype.newBuilder()
        .setVariant(ev)
        .setAlleles(Array(GenotypeAllele.Ref, GenotypeAllele.Alt).toList)
        .build()))

    assert(calls.size === 2)
  }
}
