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

import org.bdgenomics.formats.avro.{ Contig, Genotype, Variant }
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.scalatest.FunSuite

class FilterDepthSuite extends FunSuite {

  test("do not filter genotypes that have no depth info") {
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
    val filt = new DepthFilter(10)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 3)
    assert(gts.forall(_.getVariantCallingAnnotations == null))
  }

  test("do not filter genotypes that have sufficient coverage") {
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
      .setReadDepth(10)

    val seq = Seq(gt1, gt1).map(_.build())
    val filt = new DepthFilter(10)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 2)
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.contains("DEPTH>=10")))
  }

  test("filter genotypes that have low coverage") {
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
      .setReadDepth(5)

    val seq = Seq(gt1, gt1).map(_.build())
    val filt = new DepthFilter(10)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 2)
    assert(gts.forall(!_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.forall(_.getVariantCallingAnnotations.getVariantFilters.contains("DEPTH>=10")))
  }

  test("do not filter genotypes that have no depth info or that have sufficient coverage") {
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
      .setSampleId("me")
    val gt2 = Genotype.newBuilder
      .setVariant(variant)
      .setSampleId("you")
      .setReadDepth(15)

    val seq = Seq(gt1, gt1, gt2, gt2).map(_.build())
    val filt = new DepthFilter(10)
    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 4)
    assert(gts.filter(_.getVariantCallingAnnotations == null).length === 2)
    assert(gts.filter(_.getVariantCallingAnnotations != null).length === 2)
    assert(gts.filter(_.getVariantCallingAnnotations != null).forall(_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.filter(_.getVariantCallingAnnotations != null).forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.filter(_.getVariantCallingAnnotations != null).forall(_.getVariantCallingAnnotations.getVariantFilters.contains("DEPTH>=10")))
  }

  test("do not filter genotypes that have no depth info but filter low coverage calls") {
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
      .setSampleId("me")
    val gt2 = Genotype.newBuilder
      .setVariant(variant)
      .setSampleId("you")
      .setReadDepth(6)

    val seq = Seq(gt1, gt1, gt2, gt2).map(_.build())
    val filt = new DepthFilter(10)

    val gts = filt.filterGenotypes(seq)

    assert(gts.length === 4)
    assert(gts.filter(_.getVariantCallingAnnotations == null).length === 2)
    assert(gts.filter(_.getVariantCallingAnnotations == null).forall(_.getSampleId == "me"))
    assert(gts.filter(_.getVariantCallingAnnotations != null).length === 2)
    assert(gts.filter(_.getVariantCallingAnnotations != null).forall(_.getSampleId == "you"))
    assert(gts.filter(_.getVariantCallingAnnotations != null).forall(!_.getVariantCallingAnnotations.getVariantIsPassing))
    assert(gts.filter(_.getVariantCallingAnnotations != null).forall(_.getVariantCallingAnnotations.getVariantFilters.size == 1))
    assert(gts.filter(_.getVariantCallingAnnotations != null).forall(_.getVariantCallingAnnotations.getVariantFilters.contains("DEPTH>=10")))
  }
}

