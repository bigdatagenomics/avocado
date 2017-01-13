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
package org.bdgenomics.avocado.util

import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.{
  Genotype,
  GenotypeAllele,
  Variant,
  VariantCallingAnnotations
}
import scala.collection.JavaConversions._

case class TestHardFilterGenotypesArgs(
    var minQuality: Int = 30,
    var minHetSnpQualityByDepth: Float = 2.0f,
    var minHomSnpQualityByDepth: Float = 2.0f,
    var minHetIndelQualityByDepth: Float = 1.0f,
    var minHomIndelQualityByDepth: Float = 1.0f,
    var maxSnpPhredStrandBias: Float = 60.0f,
    var maxIndelPhredStrandBias: Float = 200.0f,
    var minSnpRMSMappingQuality: Float = 40.0f,
    var minIndelRMSMappingQuality: Float = 30.0f,
    var minSnpDepth: Int = 10,
    var maxSnpDepth: Int = 150,
    var minIndelDepth: Int = 20,
    var maxIndelDepth: Int = 300) extends HardFilterGenotypesArgs {
}

class HardFilterGenotypesSuite extends AvocadoFunSuite {

  def alt: Genotype.Builder = {
    Genotype.newBuilder
      .setAlleles(Seq(GenotypeAllele.ALT, GenotypeAllele.REF))
  }

  def homAlt: Genotype.Builder = {
    Genotype.newBuilder
      .setAlleles(Seq(GenotypeAllele.ALT, GenotypeAllele.ALT))
  }

  def altWithAnn: Genotype.Builder = {
    alt
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .build())
  }

  val ref = Genotype.newBuilder
    .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.REF))
    .setGenotypeQuality(99)
    .build

  test("filter out reference calls") {
    assert(!HardFilterGenotypes.filterRefCalls(ref))
    assert(HardFilterGenotypes.filterRefCalls(alt.build))
  }

  test("filter out low quality calls") {
    assert(HardFilterGenotypes.filterQuality(alt.setGenotypeQuality(50).build, 30))
    assert(!HardFilterGenotypes.filterQuality(alt.setGenotypeQuality(10).build, 30))
  }

  test("filter out genotypes for emission") {
    assert(!HardFilterGenotypes.emitGenotypeFilter(ref, 30))
    assert(HardFilterGenotypes.emitGenotypeFilter(alt.setGenotypeQuality(50).build, 30))
    assert(!HardFilterGenotypes.emitGenotypeFilter(alt.setGenotypeQuality(10).build, 30))
  }

  test("filter out genotypes with a low quality per depth") {
    val highQD = HardFilterGenotypes.hardFilterQualityByDepth(alt.setGenotypeQuality(50)
      .setReadDepth(20)
      .build, 2.0f)
    assert(highQD.isEmpty)

    val lowQD = HardFilterGenotypes.hardFilterQualityByDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 2.0f)
    assert(lowQD.isDefined)
    assert(lowQD.get === "QD<2.0")

    val homQD = HardFilterGenotypes.hardFilterQualityByDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 2.0f, hom = true)
    assert(homQD.isEmpty)

    val hetQD = HardFilterGenotypes.hardFilterQualityByDepth(homAlt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 2.0f)
    assert(homQD.isEmpty)
  }

  test("filter out genotypes with a low depth") {
    val highDP = HardFilterGenotypes.hardFilterMinDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 10)
    assert(highDP.isEmpty)

    val lowDP = HardFilterGenotypes.hardFilterMinDepth(alt.setGenotypeQuality(50)
      .setReadDepth(5)
      .build, 10)
    assert(lowDP.isDefined)
    assert(lowDP.get === "DP<10")
  }

  test("filter out genotypes with a high depth") {
    val highDP = HardFilterGenotypes.hardFilterMaxDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 10)
    assert(highDP.isDefined)
    assert(highDP.get === "DP>10")

    val lowDP = HardFilterGenotypes.hardFilterMaxDepth(alt.setGenotypeQuality(50)
      .setReadDepth(5)
      .build, 10)
    assert(lowDP.isEmpty)
  }

  test("filter out genotypes with a low RMS mapping quality") {
    val highRMQ = HardFilterGenotypes.hardFilterRMSMapQ(alt.setGenotypeQuality(50)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(50.0f)
        .build())
      .build, 30.0f)
    assert(highRMQ.isEmpty)

    val lowRMQ = HardFilterGenotypes.hardFilterRMSMapQ(alt.setGenotypeQuality(50)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(10.0f)
        .build())
      .build, 30.0f)
    assert(lowRMQ.isDefined)
    assert(lowRMQ.get === "MQ<30.0")
  }

  test("filter out genotypes with a high strand bias") {
    val highStrandBias = HardFilterGenotypes.hardFilterStrandBias(alt
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setFisherStrandBiasPValue(50.0f)
        .build())
      .build, 30.0f)
    assert(highStrandBias.isDefined)
    assert(highStrandBias.get === "FS>30.0")

    val lowStrandBias = HardFilterGenotypes.hardFilterStrandBias(alt
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setFisherStrandBiasPValue(10.0f)
        .build())
      .build, 30.0f)
    assert(lowStrandBias.isEmpty)
  }

  def validate(gt: Genotype,
               filterMsgs: Set[String] = Set.empty,
               filtersApplied: Boolean = true) {
    val optVca = Option(gt.getVariantCallingAnnotations)
    assert(optVca.isDefined)
    optVca.foreach(vca => {
      val variantIsPassing = Option(vca.getFiltersPassed)
      if (filtersApplied) {
        assert(variantIsPassing.isDefined)
        assert(variantIsPassing.get === filterMsgs.isEmpty)
        val variantMessages = vca.getFiltersFailed
        assert(variantMessages.length === filterMsgs.size)
        (0 until filterMsgs.size).foreach(idx => {
          assert(filterMsgs(variantMessages.get(idx)))
        })
      } else {
        assert(variantIsPassing.isEmpty)
      }
    })
  }

  test("update genotype where no filters were applied") {
    validate(HardFilterGenotypes.updateGenotypeWithFilters(alt.build,
      false,
      Iterable.empty), filtersApplied = false)
    validate(HardFilterGenotypes.updateGenotypeWithFilters(altWithAnn.build,
      false,
      Iterable.empty), filtersApplied = false)
  }

  test("update genotype where filters were applied and passed") {
    def validate(gt: Genotype) {
      val optVca = Option(gt.getVariantCallingAnnotations)
      assert(optVca.isDefined)
      optVca.foreach(vca => {
        val variantIsPassing = Option(vca.getFiltersPassed)
        assert(variantIsPassing.isDefined)
        assert(variantIsPassing.get)
      })
    }

    validate(HardFilterGenotypes.updateGenotypeWithFilters(alt.build,
      true,
      Iterable.empty))
    validate(HardFilterGenotypes.updateGenotypeWithFilters(altWithAnn.build,
      true,
      Iterable.empty))
  }

  test("update genotype where filters were applied and failed") {
    validate(HardFilterGenotypes.updateGenotypeWithFilters(alt.build,
      true,
      Iterable("NOT_VARIANT")),
      filterMsgs = Set("NOT_VARIANT"))
    validate(HardFilterGenotypes.updateGenotypeWithFilters(altWithAnn.build,
      true,
      Iterable("SUCH_REF", "NOT_VARIANT")),
      filterMsgs = Set("SUCH_REF", "NOT_VARIANT"))
  }

  val defaultSnpFilters = HardFilterGenotypes.buildSnpHardFilters(
    TestHardFilterGenotypesArgs())
  val defaultIndelFilters = HardFilterGenotypes.buildIndelHardFilters(
    TestHardFilterGenotypesArgs())

  test("discard a ref genotype call") {
    val optGt = HardFilterGenotypes.filterGenotype(ref,
      20,
      Iterable.empty,
      Iterable.empty)
    assert(optGt.isEmpty)
  }

  test("discard a genotype whose quality is too low") {
    val lowQualAlt = alt.setGenotypeQuality(10).build
    val optGt = HardFilterGenotypes.filterGenotype(lowQualAlt,
      20,
      Iterable.empty,
      Iterable.empty)
    assert(optGt.isEmpty)
  }

  test("build filters and apply to snp") {
    val failingSnp = alt.setVariant(Variant.newBuilder
      .setReferenceAllele("A")
      .setAlternateAllele("T")
      .build)
      .setGenotypeQuality(30)
      .setReadDepth(20)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(35.0f)
        .setFisherStrandBiasPValue(70.0f)
        .build())
      .build()

    val optGt = HardFilterGenotypes.filterGenotype(failingSnp,
      20,
      defaultSnpFilters,
      defaultIndelFilters)
    assert(optGt.isDefined)
    validate(optGt.get,
      filterMsgs = Set("MQ<40.0", "QD<2.0", "FS>60.0"))
  }

  test("build filters and apply to indel") {
    val passingIndel = alt.setVariant(Variant.newBuilder
      .setReferenceAllele("A")
      .setAlternateAllele("AT")
      .build)
      .setGenotypeQuality(30)
      .setReadDepth(20)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(35.0f)
        .setFisherStrandBiasPValue(70.0f)
        .build())
      .build()

    val optPassingGt = HardFilterGenotypes.filterGenotype(passingIndel,
      20,
      defaultSnpFilters,
      defaultIndelFilters)
    assert(optPassingGt.isDefined)
    validate(optPassingGt.get)

    val failingIndel = alt.setVariant(Variant.newBuilder
      .setReferenceAllele("A")
      .setAlternateAllele("AT")
      .build)
      .setGenotypeQuality(15)
      .setReadDepth(20)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(25.0f)
        .setFisherStrandBiasPValue(250.0f)
        .build())
      .build()

    val optFailingGt = HardFilterGenotypes.filterGenotype(failingIndel,
      10,
      defaultSnpFilters,
      defaultIndelFilters)
    assert(optFailingGt.isDefined)
    validate(optFailingGt.get,
      filterMsgs = Set("MQ<30.0", "QD<1.0", "FS>200.0"))
  }
}
