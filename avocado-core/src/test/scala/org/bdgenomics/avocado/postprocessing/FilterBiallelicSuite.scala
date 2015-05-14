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

import org.bdgenomics.formats.avro.{ Contig, Genotype, Variant, VariantCallingAnnotations }
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.scalatest.FunSuite

class FilterBiallelicSuite extends FunSuite {

  test("genotypes that do not have the multiallelic field set should pass") {
    val contig = Contig.newBuilder
      .setContigName("0")
      .build
    val variant = Variant.newBuilder
      .setContig(contig)
      .setStart(0)
      .setEnd(1)
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .build
    val gt1 = Genotype.newBuilder
      .setVariant(variant)

    val seq = Seq(gt1, gt1, gt1).map(_.build())
    val filt = new BiallelicFilter()
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 3)
    assert(gts.forall(_.getVariantCallingAnnotations != null))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.contains("BIALLELIC")))
  }

  test("do not filter genotypes that are biallelic") {
    val contig = Contig.newBuilder
      .setContigName("0")
      .build
    val variant = Variant.newBuilder
      .setContig(contig)
      .setStart(0)
      .setEnd(1)
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .build
    val gt1 = Genotype.newBuilder
      .setVariant(variant)
      .setSplitFromMultiAllelic(false)

    val seq = Seq(gt1, gt1).map(_.build())
    val filt = new BiallelicFilter()
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 2)
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.contains("BIALLELIC")))
  }

  test("filter genotypes that are from a multiallelic site") {
    val contig = Contig.newBuilder
      .setContigName("0")
      .build
    val variant = Variant.newBuilder
      .setContig(contig)
      .setStart(0)
      .setEnd(1)
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .build
    val gt1 = Genotype.newBuilder
      .setVariant(variant)
      .setSplitFromMultiAllelic(true)

    val seq = Seq(gt1, gt1).map(_.build())
    val filt = new BiallelicFilter()
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 2)
    assert(gts.forall(!_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.contains("BIALLELIC")))
  }
}

