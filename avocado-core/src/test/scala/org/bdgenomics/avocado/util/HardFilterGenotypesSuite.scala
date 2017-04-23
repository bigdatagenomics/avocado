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

import org.bdgenomics.adam.rdd.ADAMContext._
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
    var maxIndelDepth: Int = 300,
    var minHetSnpAltAllelicFraction: Float = 0.333f,
    var maxHetSnpAltAllelicFraction: Float = 0.666f,
    var minHomSnpAltAllelicFraction: Float = 0.666f,
    var minHetIndelAltAllelicFraction: Float = 0.333f,
    var maxHetIndelAltAllelicFraction: Float = 0.666f,
    var minHomIndelAltAllelicFraction: Float = 0.666f) extends HardFilterGenotypesArgs {
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
    .setVariant(Variant.newBuilder
      .setReferenceAllele("A")
      .setAlternateAllele("T")
      .build)
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
    assert(!HardFilterGenotypes.emitGenotypeFilter(ref, 30, true))
    assert(HardFilterGenotypes.emitGenotypeFilter(ref, 30, false))
    assert(HardFilterGenotypes.emitGenotypeFilter(alt.setGenotypeQuality(50).build, 30, true))
    assert(!HardFilterGenotypes.emitGenotypeFilter(alt.setGenotypeQuality(10).build, 30, true))
  }

  test("filter out genotypes with a low quality per depth") {
    val highQD = HardFilterGenotypes.hardFilterQualityByDepth(alt.setGenotypeQuality(50)
      .setReadDepth(20)
      .build, 2.0f, "QD")
    assert(highQD.isEmpty)

    val lowQD = HardFilterGenotypes.hardFilterQualityByDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 2.0f, "QD")
    assert(lowQD.isDefined)
    assert(lowQD.get === "QD")

    val homQD = HardFilterGenotypes.hardFilterQualityByDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 2.0f, "QD", hom = true)
    assert(homQD.isEmpty)

    val hetQD = HardFilterGenotypes.hardFilterQualityByDepth(homAlt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 2.0f, "QD")
    assert(homQD.isEmpty)
  }

  test("filter out genotypes with a low depth") {
    val highDP = HardFilterGenotypes.hardFilterMinDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 10, "MINDP")
    assert(highDP.isEmpty)

    val lowDP = HardFilterGenotypes.hardFilterMinDepth(alt.setGenotypeQuality(50)
      .setReadDepth(5)
      .build, 10, "MINDP")
    assert(lowDP.isDefined)
    assert(lowDP.get === "MINDP")
  }

  test("filter out genotypes with a high depth") {
    val highDP = HardFilterGenotypes.hardFilterMaxDepth(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .build, 10, "MAXDP")
    assert(highDP.isDefined)
    assert(highDP.get === "MAXDP")

    val lowDP = HardFilterGenotypes.hardFilterMaxDepth(alt.setGenotypeQuality(50)
      .setReadDepth(5)
      .build, 10, "MAXDP")
    assert(lowDP.isEmpty)
  }

  test("filter out genotypes with a low RMS mapping quality") {
    val highRMQ = HardFilterGenotypes.hardFilterRMSMapQ(alt.setGenotypeQuality(50)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(50.0f)
        .build())
      .build, 30.0f, "MQ")
    assert(highRMQ.isEmpty)

    val lowRMQ = HardFilterGenotypes.hardFilterRMSMapQ(alt.setGenotypeQuality(50)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(10.0f)
        .build())
      .build, 30.0f, "MQ")
    assert(lowRMQ.isDefined)
    assert(lowRMQ.get === "MQ")
  }

  test("filter out genotypes with a high strand bias") {
    val highStrandBias = HardFilterGenotypes.hardFilterStrandBias(alt
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setFisherStrandBiasPValue(50.0f)
        .build())
      .build, 30.0f, "FS")
    assert(highStrandBias.isDefined)
    assert(highStrandBias.get === "FS")

    val lowStrandBias = HardFilterGenotypes.hardFilterStrandBias(alt
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setFisherStrandBiasPValue(10.0f)
        .build())
      .build, 30.0f, "FS")
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
      Iterable.empty,
      true)
    assert(optGt.isEmpty)
  }

  test("keep a ref genotype call") {
    val optGt = HardFilterGenotypes.filterGenotype(ref,
      20,
      Iterable.empty,
      Iterable.empty,
      false)
    assert(optGt.isDefined)
  }

  test("discard a genotype whose quality is too low") {
    val lowQualAlt = alt.setGenotypeQuality(10).build
    val optGt = HardFilterGenotypes.filterGenotype(lowQualAlt,
      20,
      Iterable.empty,
      Iterable.empty,
      true)
    assert(optGt.isEmpty)
  }

  test("build filters and apply to snp") {
    val failingSnp = alt.setVariant(Variant.newBuilder
      .setReferenceAllele("A")
      .setAlternateAllele("T")
      .build)
      .setGenotypeQuality(30)
      .setAlternateReadDepth(15)
      .setReadDepth(20)
      .setVariantCallingAnnotations(VariantCallingAnnotations.newBuilder
        .setRmsMapQ(35.0f)
        .setFisherStrandBiasPValue(70.0f)
        .build())
      .build()

    val optGt = HardFilterGenotypes.filterGenotype(failingSnp,
      20,
      defaultSnpFilters,
      defaultIndelFilters,
      true)
    assert(optGt.isDefined)
    validate(optGt.get,
      filterMsgs = Set("SNPMQ", "HETSNPQD", "SNPFS", "HETSNPMAXAF"))
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
      defaultIndelFilters,
      true)
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
      defaultIndelFilters,
      true)
    assert(optFailingGt.isDefined)
    validate(optFailingGt.get,
      filterMsgs = Set("INDELMQ", "HETINDELQD", "INDELFS"))
  }

  sparkTest("test adding filters") {
    val path = testFile("random.vcf")
    val initialRdd = sc.loadGenotypes(path)
    val filteredRdd = HardFilterGenotypes(
      initialRdd,
      TestHardFilterGenotypesArgs())
    assert(filteredRdd.headerLines.size === initialRdd.headerLines.size + 18)
  }

  test("filter out genotypes with a low allelic fraction") {
    val highAF = HardFilterGenotypes.hardFilterMinAllelicFraction(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .setAlternateReadDepth(25)
      .build, 0.333f, "MINAF", hom = false)
    assert(highAF.isEmpty)

    val lowAF = HardFilterGenotypes.hardFilterMinAllelicFraction(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .setAlternateReadDepth(10)
      .build, 0.333f, "MINAF", hom = false)
    assert(lowAF.isDefined)
    assert(lowAF.get === "MINAF")

    val lowAFHet = HardFilterGenotypes.hardFilterMinAllelicFraction(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .setAlternateReadDepth(10)
      .build, 0.333f, "MINAF", hom = true)
    assert(lowAFHet.isEmpty)
  }

  test("filter out genotypes with a high allelic fraction") {
    val highAF = HardFilterGenotypes.hardFilterMaxAllelicFraction(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .setAlternateReadDepth(45)
      .build, 0.666f, "MAXAF")
    assert(highAF.isDefined)
    assert(highAF.get === "MAXAF")

    val lowAF = HardFilterGenotypes.hardFilterMaxAllelicFraction(alt.setGenotypeQuality(50)
      .setReadDepth(50)
      .setAlternateReadDepth(30)
      .build, 0.666f, "MAXAF")
    assert(lowAF.isEmpty)

    val hom = HardFilterGenotypes.hardFilterMaxAllelicFraction(homAlt.setGenotypeQuality(50)
      .setReadDepth(50)
      .setAlternateReadDepth(45)
      .build, 0.666f, "MAXAF")
    assert(hom.isEmpty)
  }
}
