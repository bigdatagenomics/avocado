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

class FilterStrandBiasSuite extends FunSuite {

  test("do not filter genotypes that have no variant calling annotations") {
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
    val filt = new StrandBiasFilter(10.0f)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 3)
    assert(gts.forall(_.getVariantCallingAnnotations == null))
  }

  test("do not filter genotypes that have no strand bias info") {
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

    val seq = Seq(gt1, gt1, gt1).map(v => {
      v.setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder()
        .build()).build()
    })
    val filt = new StrandBiasFilter(10.0f)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 3)
    assert(gts.forall(_.getVariantCallingAnnotations != null))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 0))
  }

  test("do not filter genotypes that have strand bias phred score above threshold") {
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

    val seq = Seq(gt1, gt1).map(v => {
      v.setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder()
        .setFisherStrandBiasPValue(11.0f)
        .build()).build()
    })
    val filt = new StrandBiasFilter(10.0f)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 2)
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.contains("FISHER_STRAND_BIAS_PVALUE>=10.000000")))
  }

  test("filter genotypes that have strand bias phred score below threshold") {
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

    val seq = Seq(gt1, gt1).map(v => {
      v.setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder()
        .setFisherStrandBiasPValue(9.0f)
        .build()).build()
    })
    val filt = new StrandBiasFilter(10.0f)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 2)
    assert(gts.forall(!_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.contains("FISHER_STRAND_BIAS_PVALUE>=10.000000")))
  }
}

