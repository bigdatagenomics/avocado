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

import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.{
  Variant,
  Genotype
}
import org.bdgenomics.utils.misc.MathUtils
import scala.collection.JavaConversions._

class SquareOffReferenceModelSuite extends AvocadoFunSuite {

  test("don't trim a snp") {
    assert(SquareOffReferenceModel.trimRight("A",
      "T") === 0)
  }

  test("trim a mnp") {
    assert(SquareOffReferenceModel.trimRight("ACG",
      "TAG") === 1)
  }

  test("trim an insertion") {
    assert(SquareOffReferenceModel.trimRight("TAA",
      "TCAA") === 2)
  }

  test("don't trim a deletion") {
    assert(SquareOffReferenceModel.trimRight("TCCAGT",
      "TG") === 0)
  }

  sparkTest("extract variants finds sites with a called alt") {
    val gvcf = resourceUrl("gvcf_multiallelic.g.vcf")

    val allVariants = sc.loadGenotypes(gvcf.toString)

    val discoveredVariants = SquareOffReferenceModel.extractVariants(
      allVariants)

    val variants = discoveredVariants.rdd.collect
    assert(variants.size === 3)
    assert(variants.forall(_.getReferenceName == "chr22"))
    val s602 = variants.filter(_.getStart == 16157602L)
    assert(s602.size === 1)
    assert(s602.forall(_.getReferenceAllele == "G"))
    assert(s602.forall(_.getAlternateAllele == "C"))
    assert(s602.forall(_.getEnd == 16157603L))
    val s095 = variants.filter(_.getStart == 18030095L)
    assert(s095.size === 2)
    assert(s095.forall(_.getAlternateAllele == "T"))
    val s095taa = variants.filter(_.getReferenceAllele == "TAA")
    assert(s095taa.size === 1)
    assert(s095taa.forall(_.getEnd == 18030098L))
    val s095ta = variants.filter(_.getReferenceAllele == "TA")
    assert(s095ta.size == 1)
    assert(s095ta.forall(_.getEnd == 18030097L))
  }

  test("find genotype if variant is present") {
    val variant = Variant.newBuilder()
      .setReferenceName("ctg")
      .setStart(1L)
      .setEnd(2L)
      .setReferenceAllele("A")
      .setAlternateAllele("T")
      .build
    val genotypes = Iterable(Genotype.newBuilder
      .setReferenceName("ctg")
      .setStart(1L)
      .setEnd(2L)
      .setVariant(Variant.newBuilder()
        .setReferenceAllele("A")
        .setAlternateAllele("T")
        .build)
      .build)

    val optFoundGenotype = SquareOffReferenceModel.hasScoredVariant(variant,
      genotypes)

    assert(optFoundGenotype === genotypes.headOption)
  }

  test("don't find genotype if variant is not present") {
    val variant = Variant.newBuilder()
      .setReferenceName("ctg")
      .setStart(1L)
      .setEnd(2L)
      .setReferenceAllele("A")
      .setAlternateAllele("T")
      .build
    val genotypes = Iterable(Genotype.newBuilder
      .setReferenceName("ctg")
      .setStart(1L)
      .setEnd(10L)
      .setVariant(Variant.newBuilder()
        .setReferenceAllele("A")
        .build)
      .build)

    val optFoundGenotype = SquareOffReferenceModel.hasScoredVariant(variant,
      genotypes)

    assert(optFoundGenotype.isEmpty)
  }

  test("excise a genotype from a reference block") {
    val variant = Variant.newBuilder
      .setStart(100L)
      .setEnd(101L)
      .setReferenceName("ctg")
      .setReferenceAllele("A")
      .setAlternateAllele("G")
      .build
    val genotypes = Iterable(Genotype.newBuilder
      .setStart(90L)
      .setEnd(110L)
      .setReferenceName("ctg")
      .setNonReferenceLikelihoods(Seq(0.0, -1.0, -2.0)
        .map(d => d: java.lang.Double))
      .build)

    val optExcisedGenotype = SquareOffReferenceModel.exciseGenotypeFromReferenceModel(
      variant,
      genotypes)

    assert(optExcisedGenotype.isDefined)
    optExcisedGenotype.foreach(gt => {
      assert(gt.getStart === 100L)
      assert(gt.getEnd === 101L)
      assert(gt.getReferenceName === "ctg")
      assert(gt.getVariant.getReferenceAllele === "A")
      assert(gt.getVariant.getAlternateAllele === "G")
      assert(gt.getGenotypeLikelihoods.size === 3)
      assert(MathUtils.fpEquals(gt.getGenotypeLikelihoods.get(0), 0.0))
      assert(MathUtils.fpEquals(gt.getGenotypeLikelihoods.get(1), -1.0))
      assert(MathUtils.fpEquals(gt.getGenotypeLikelihoods.get(2), -2.0))
    })
  }

  test("square off a site with data from multiple samples") {
    val variant = Variant.newBuilder
      .setStart(100L)
      .setEnd(101L)
      .setReferenceName("ctg")
      .setReferenceAllele("A")
      .setAlternateAllele("G")
      .build
    val genotypes = Iterable(Genotype.newBuilder
      .setReferenceName("ctg")
      .setStart(100L)
      .setEnd(101L)
      .setVariant(Variant.newBuilder()
        .setReferenceAllele("A")
        .setAlternateAllele("G")
        .build)
      .setSampleId("sample1")
      .build, Genotype.newBuilder
      .setStart(90L)
      .setEnd(110L)
      .setReferenceName("ctg")
      .setNonReferenceLikelihoods(Seq(0.0, -1.0, -2.0)
        .map(d => d: java.lang.Double))
      .setSampleId("sample2")
      .build)

    val vc = SquareOffReferenceModel.squareOffSite(variant,
      genotypes)

    assert(vc.variant.variant === variant)
    assert(vc.genotypes.size === 2)
    assert(vc.genotypes.forall(_.getStart == 100L))
    assert(vc.genotypes.forall(_.getEnd == 101L))
    assert(vc.genotypes.forall(_.getReferenceName == "ctg"))
    assert(vc.genotypes.forall(_.getVariant.getReferenceAllele == "A"))
    assert(vc.genotypes.forall(_.getVariant.getAlternateAllele == "G"))
    assert(vc.genotypes.count(_.getSampleId == "sample1") === 1)
    val s2gts = vc.genotypes.filter(_.getSampleId == "sample2")
    assert(s2gts.size === 1)
    val s2gt = s2gts.head
    assert(MathUtils.fpEquals(s2gt.getGenotypeLikelihoods.get(0), 0.0))
    assert(MathUtils.fpEquals(s2gt.getGenotypeLikelihoods.get(1), -1.0))
    assert(MathUtils.fpEquals(s2gt.getGenotypeLikelihoods.get(2), -2.0))
  }
}
