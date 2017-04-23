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

import htsjdk.samtools.ValidationStringency
import org.bdgenomics.adam.models.{
  RecordGroup,
  RecordGroupDictionary,
  SequenceDictionary,
  VariantContext
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.{
  Genotype,
  GenotypeAllele,
  Variant
}
import scala.collection.JavaConversions._

class TrioCallerSuite extends AvocadoFunSuite {

  def makeRdd(recordGroups: RecordGroupDictionary): AlignmentRecordRDD = {
    AlignmentRecordRDD(sc.emptyRDD,
      SequenceDictionary.empty,
      recordGroups)
  }

  sparkTest("cannot have a sample with no record groups") {
    intercept[IllegalArgumentException] {
      TrioCaller.extractSampleId(makeRdd(RecordGroupDictionary.empty))
    }
  }

  sparkTest("cannot have a sample with discordant sample ids") {
    intercept[IllegalArgumentException] {
      TrioCaller.extractSampleId(makeRdd(RecordGroupDictionary(Seq(
        RecordGroup("sample1", "rg1"),
        RecordGroup("sample2", "rg2")))))
    }
  }

  sparkTest("extract id from a single read group") {
    val sampleId = TrioCaller.extractSampleId(makeRdd(RecordGroupDictionary(Seq(
      RecordGroup("sample1", "rg1")))))
    assert(sampleId === "sample1")
  }

  sparkTest("extract id from multiple read groups") {
    val sampleId = TrioCaller.extractSampleId(makeRdd(RecordGroupDictionary(Seq(
      RecordGroup("sample1", "rg1"),
      RecordGroup("sample1", "rg2")))))
    assert(sampleId === "sample1")
  }

  val variant = Variant.newBuilder
    .setContigName("chr")
    .setStart(100L)
    .setEnd(101L)
    .setReferenceAllele("A")
    .setAlternateAllele("T")
    .build

  test("filter an empty site") {
    val emptyVc = VariantContext(variant, Iterable.empty)
    assert(TrioCaller.filterRef(emptyVc))
  }

  test("filter a site with only ref calls") {
    val refVc = VariantContext(variant, Iterable(Genotype.newBuilder
      .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.REF))
      .build()))
    assert(TrioCaller.filterRef(refVc))
  }

  test("keep a site with a non-ref call") {
    val nonRefVc = VariantContext(variant, Iterable(Genotype.newBuilder
      .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.REF))
      .build(),
      Genotype.newBuilder
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .build))
    assert(!TrioCaller.filterRef(nonRefVc))
  }

  test("fill in no-calls for site with missing parents") {
    val emptyVc = VariantContext(variant, Iterable(Genotype.newBuilder
      .setSampleId("proband")
      .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
      .build))
    val vc = TrioCaller.processVariant(emptyVc,
      "parent1",
      "parent2",
      "proband")
    assert(vc.genotypes.size === 3)
    vc.genotypes.foreach(gt => {
      assert(gt.getAlleles.size === 2)
      if (gt.getSampleId == "proband") {
        assert(gt.getAlleles.count(_ == GenotypeAllele.REF) === 1)
        assert(gt.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
      } else {
        assert(gt.getAlleles.forall(_ == GenotypeAllele.NO_CALL))
      }
    })
    assert(vc.genotypes.count(_.getSampleId == "parent1") === 1)
    assert(vc.genotypes.count(_.getSampleId == "parent2") === 1)
    assert(vc.genotypes.count(_.getSampleId == "proband") === 1)
  }

  test("pass through site with odd copy number") {
    val mixedCnVc = VariantContext(variant, Iterable(Genotype.newBuilder
      .setSampleId("proband")
      .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
      .build,
      Genotype.newBuilder
        .setSampleId("parent1")
        .setAlleles(Seq(GenotypeAllele.REF))
        .build,
      Genotype.newBuilder
        .setSampleId("parent2")
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .build))
    val vc = TrioCaller.processVariant(mixedCnVc,
      "parent1",
      "parent2",
      "proband")
    assert(vc.variant.variant === mixedCnVc.variant.variant)
    val gtSet = mixedCnVc.genotypes.toSet
    vc.genotypes.foreach(gt => assert(gtSet(gt)))
  }

  test("confirm call at site where proband and parents are consistent and phase") {
    val conVc = VariantContext(variant, Iterable(Genotype.newBuilder
      .setSampleId("proband")
      .setAlleles(Seq(GenotypeAllele.ALT, GenotypeAllele.REF))
      .build,
      Genotype.newBuilder
        .setSampleId("parent1")
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.REF))
        .build,
      Genotype.newBuilder
        .setSampleId("parent2")
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .build))
    val vc = TrioCaller.processVariant(conVc,
      "parent1",
      "parent2",
      "proband")
    assert(vc.genotypes.size === 3)
    assert(vc.genotypes.count(_.getSampleId == "parent1") === 1)
    assert(vc.genotypes.count(_.getSampleId == "parent2") === 1)
    assert(vc.genotypes.count(_.getSampleId == "proband") === 1)
    vc.genotypes.filter(_.getSampleId == "proband")
      .foreach(gt => {
        assert(gt.getAlleles.size === 2)
        assert(gt.getAlleles.get(0) === GenotypeAllele.REF)
        assert(gt.getAlleles.get(1) === GenotypeAllele.ALT)
        assert(gt.getPhased)
      })
  }

  test("confirm call at site where proband and parents are consistent but cannot phase") {
    val conVc = VariantContext(variant, Iterable(Genotype.newBuilder
      .setSampleId("proband")
      .setAlleles(Seq(GenotypeAllele.ALT, GenotypeAllele.REF))
      .build,
      Genotype.newBuilder
        .setSampleId("parent1")
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .build,
      Genotype.newBuilder
        .setSampleId("parent2")
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .build))
    val vc = TrioCaller.processVariant(conVc,
      "parent1",
      "parent2",
      "proband")
    assert(vc.genotypes.size === 3)
    assert(vc.genotypes.count(_.getSampleId == "parent1") === 1)
    assert(vc.genotypes.count(_.getSampleId == "parent2") === 1)
    assert(vc.genotypes.count(_.getSampleId == "proband") === 1)
    vc.genotypes.filter(_.getSampleId == "proband")
      .foreach(gt => {
        assert(gt.getAlleles.size === 2)
        assert(gt.getAlleles.count(_ == GenotypeAllele.REF) === 1)
        assert(gt.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
        assert(!gt.getPhased)
      })
  }

  test("invalidate call at site where proband and parents are inconsistent") {
    val incVc = VariantContext(variant, Iterable(Genotype.newBuilder
      .setSampleId("proband")
      .setAlleles(Seq(GenotypeAllele.ALT, GenotypeAllele.ALT))
      .build,
      Genotype.newBuilder
        .setSampleId("parent1")
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.REF))
        .build,
      Genotype.newBuilder
        .setSampleId("parent2")
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .build))
    val vc = TrioCaller.processVariant(incVc,
      "parent1",
      "parent2",
      "proband")
    assert(vc.genotypes.size === 3)
    assert(vc.genotypes.count(_.getSampleId == "parent1") === 1)
    assert(vc.genotypes.count(_.getSampleId == "parent2") === 1)
    assert(vc.genotypes.count(_.getSampleId == "proband") === 1)
    vc.genotypes.filter(_.getSampleId == "proband")
      .foreach(gt => {
        assert(gt.getAlleles.size === 2)
        assert(gt.getAlleles.count(_ == GenotypeAllele.NO_CALL) == 2)
      })
  }

  sparkTest("end-to-end trio call test") {
    val inputVcf = resourceUrl("trio.1_837214.vcf")
    val compVcf = resourceUrl("trio.1_837214.phased.vcf")
    val outputVcf = tmpFile("trio.vcf")
    println(outputVcf)

    val gts = sc.loadGenotypes(inputVcf.toString)

    TrioCaller(gts, "NA12891", "NA12892", "NA12878")
      .saveAsVcf(outputVcf, true, ValidationStringency.STRICT)

    checkFiles(outputVcf, compVcf.getPath)
  }
}
