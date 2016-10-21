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
import org.bdgenomics.avocado.models.Observation
import org.bdgenomics.avocado.util.{
  HardFilterGenotypes,
  TestHardFilterGenotypesArgs
}
import org.bdgenomics.formats.avro.{
  AlignmentRecord,
  GenotypeAllele,
  Variant
}
import org.bdgenomics.utils.misc.MathUtils
import scala.collection.JavaConversions._
import scala.math.log

class BiallelicGenotyperSuite extends AvocadoFunSuite {

  val os = new ObserverSuite()

  val perfectRead = AlignmentRecord.newBuilder
    .setContigName("1")
    .setStart(10L)
    .setEnd(25L)
    .setCigar("15M")
    .setMismatchingPositions("15")
    .setSequence("ATGGTCCACGAATAA")
    .setQual("DEFGHIIIIIHGFED")
    .setMapq(50)
    .setReadMapped(true)
    .build

  val snpRead = AlignmentRecord.newBuilder(perfectRead)
    .setMismatchingPositions("6A8")
    .setSequence("ATGGTCAACGAATAA")
    .setMapq(40)
    .setReadNegativeStrand(true)
    .build

  val snp = Variant.newBuilder
    .setContigName("1")
    .setStart(16L)
    .setEnd(17L)
    .setReferenceAllele("C")
    .setAlternateAllele("A")
    .build

  test("scoring read that overlaps no variants should return empty observations") {
    val perfectScores = BiallelicGenotyper.readToObservations((perfectRead, Iterable.empty), 2)
    assert(perfectScores.isEmpty)

    val snpScores = BiallelicGenotyper.readToObservations((snpRead, Iterable.empty), 2)
    assert(snpScores.isEmpty)
  }

  test("score snp in a read with no evidence of the snp") {
    val scores = BiallelicGenotyper.readToObservations(
      (perfectRead, Iterable(snp)), 2)
    assert(scores.size === 1)

    val (snpVariant, snpObservation) = scores.head
    assert(snpVariant === snp)
    assert(snpObservation.squareMapQ === 50 * 50)
    assert(snpObservation.alleleCoverage === 0)
    assert(snpObservation.otherCoverage === 1)
    assert(snpObservation.alleleForwardStrand === 0)
    assert(snpObservation.otherForwardStrand === 1)
    assert(snpObservation.copyNumber === 2)
    assert(MathUtils.fpEquals(snpObservation.otherLogLikelihoods(0), os.logL(0, 2, 0.9999, 0.99999)))
    assert(MathUtils.fpEquals(snpObservation.otherLogLikelihoods(1), os.logL(1, 2, 0.9999, 0.99999)))
    assert(MathUtils.fpEquals(snpObservation.otherLogLikelihoods(2), os.logL(2, 2, 0.9999, 0.99999)))
  }

  test("score snp in a read with evidence of the snp") {
    val scores = BiallelicGenotyper.readToObservations(
      (snpRead, Iterable(snp)), 2)
    assert(scores.size === 1)

    val (snpVariant, snpObservation) = scores.head
    assert(snpVariant === snp)
    assert(snpObservation.squareMapQ === 40 * 40)
    assert(snpObservation.alleleCoverage === 1)
    assert(snpObservation.otherCoverage === 0)
    assert(snpObservation.alleleForwardStrand === 0)
    assert(snpObservation.otherForwardStrand === 0)
    assert(snpObservation.copyNumber === 2)
    assert(MathUtils.fpEquals(snpObservation.alleleLogLikelihoods(0), os.logL(0, 2, 0.9999, 0.9999)))
    assert(MathUtils.fpEquals(snpObservation.alleleLogLikelihoods(1), os.logL(1, 2, 0.9999, 0.9999)))
    assert(MathUtils.fpEquals(snpObservation.alleleLogLikelihoods(2), os.logL(2, 2, 0.9999, 0.9999)))
  }

  test("build genotype for het snp") {
    val obs = Observation(3, 5,
      40 * 40 * 16,
      Array(-10.0, -1.0, -10.0),
      Array(-12.0, -4.0, -24.0),
      7, 9)
    val genotype = BiallelicGenotyper.observationToGenotype((snp, obs), "sample")

    assert(genotype.getVariant === snp)
    assert(genotype.getStart === snp.getStart)
    assert(genotype.getEnd === snp.getEnd)
    assert(genotype.getContigName === snp.getContigName)
    assert(genotype.getSampleId === "sample")
    assert(genotype.getGenotypeQuality === 36)
    assert(genotype.getAlleles.size === 2)
    assert(genotype.getAlleles.get(0) === GenotypeAllele.Alt)
    assert(genotype.getAlleles.get(1) === GenotypeAllele.Ref)
    assert(genotype.getReferenceReadDepth === 9)
    assert(genotype.getAlternateReadDepth === 7)
    assert(genotype.getStrandBiasComponents.size === 4)
    assert(genotype.getStrandBiasComponents.get(0) === 5)
    assert(genotype.getStrandBiasComponents.get(1) === 4)
    assert(genotype.getStrandBiasComponents.get(2) === 3)
    assert(genotype.getStrandBiasComponents.get(3) === 4)
    assert(genotype.getGenotypeLikelihoods.size === 3)
  }

  sparkTest("force call possible STR/indel") {
    val readPath = resourceUrl("NA12878.chr1.104160.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val variants = DiscoverVariants(reads)
      .transform(rdd => {
        rdd.filter(_.getAlternateAllele.startsWith("AACAC"))
      })
    assert(variants.rdd.count === 3)

    val genotypes = BiallelicGenotyper.call(reads,
      variants,
      2,
      optDesiredPartitionCount = Some(26))
    val gts = genotypes.rdd.collect
    assert(gts.size === 3)
    assert(gts.map(_.getReadDepth).toSet.size === 1)

    assert(gts.count(gt => {
      !gt.getAlleles.forall(_ == GenotypeAllele.Ref)
    }) === 2)
    assert(gts.sortBy(gt => gt.getGenotypeQuality)
      .reverse // sortBy doesn't have ascending/descending options!
      .head
      .getVariant
      .getAlternateAllele
      .length === 9)
  }

  test("log space factorial") {
    val fact0 = BiallelicGenotyper.logFactorial(0)
    assert(MathUtils.fpEquals(fact0, 0.0))

    val fact1 = BiallelicGenotyper.logFactorial(1)
    assert(MathUtils.fpEquals(fact1, 0.0))

    val fact5 = BiallelicGenotyper.logFactorial(5)
    assert(MathUtils.fpEquals(fact5, log(120.0)))

    // overflow at 13!
    val fact10 = BiallelicGenotyper.logFactorial(13)
    assert(MathUtils.fpEquals(fact10, log(6227020800L.toDouble)))
  }

  test("fisher test for strand bias") {
    val phredPBias1 = BiallelicGenotyper.fisher(43, 17, 2, 7)

    assert(MathUtils.fpEquals(phredPBias1, 22.185272))

    val phredPBias2 = BiallelicGenotyper.fisher(0, 12, 10, 2)
    assert(MathUtils.fpEquals(phredPBias2, 44.729904))
  }

  sparkTest("discover and call simple SNP") {
    val readPath = resourceUrl("NA12878_snp_A2G_chr20_225058.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val genotypes = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optDesiredPartitionCount = Some(26))
      .transform(rdd => {
        rdd.filter(gt => {
          gt.getReadDepth > 10 &&
            !gt.getAlleles.forall(_ == GenotypeAllele.Ref) &&
            gt.getGenotypeQuality > 30
        })
      })
    val gts = genotypes.rdd.collect
    assert(gts.size === 1)
    val gt = gts.head
    assert(gt.getVariant.getStart === 225057L)
    assert(gt.getVariant.getEnd === 225058L)
    assert(gt.getVariant.getReferenceAllele === "A")
    assert(gt.getVariant.getAlternateAllele === "G")
    assert(gt.getAlleles.count(_ == GenotypeAllele.Ref) === 1)
    assert(gt.getAlleles.count(_ == GenotypeAllele.Alt) === 1)
  }

  sparkTest("discover and call short indel") {
    val readPath = resourceUrl("NA12878.chr1.832736.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val genotypes = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optDesiredPartitionCount = Some(26))
    val filteredGenotypes = HardFilterGenotypes(genotypes,
      TestHardFilterGenotypesArgs())
      .transform(rdd => {
        rdd.filter(_.getVariantCallingAnnotations
          .getVariantIsPassing)
      })

    val gts = filteredGenotypes.rdd.collect

    assert(gts.size === 1)
    val gt = gts.head
    assert(gt.getVariant.getStart === 832735L)
    assert(gt.getVariant.getEnd === 832736L)
    assert(gt.getVariant.getReferenceAllele === "A")
    assert(gt.getVariant.getAlternateAllele === "AGTTTT")
    assert(gt.getAlleles.count(_ == GenotypeAllele.Ref) === 1)
    assert(gt.getAlleles.count(_ == GenotypeAllele.Alt) === 1)
  }

  sparkTest("discover and call het and hom snps") {
    val readPath = resourceUrl("NA12878.chr1.839395.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val genotypes = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optDesiredPartitionCount = Some(26))
    val filteredGenotypes = HardFilterGenotypes(genotypes,
      TestHardFilterGenotypesArgs())
      .transform(rdd => {
        rdd.filter(gt => {
          gt.getVariantCallingAnnotations
            .getVariantIsPassing &&
            gt.getReadDepth > 30
        })
      })

    val gts = filteredGenotypes.rdd.collect
    assert(gts.size === 2)

    val hetGt = gts.filter(gt => {
      gt.getAlleles.count(_ == GenotypeAllele.Alt) == 1
    }).head
    assert(hetGt.getVariant.getStart === 839404L)
    assert(hetGt.getVariant.getEnd === 839405L)
    assert(hetGt.getVariant.getReferenceAllele === "G")
    assert(hetGt.getVariant.getAlternateAllele === "A")

    val homGt = gts.filter(gt => {
      gt.getAlleles.count(_ == GenotypeAllele.Alt) == 2
    }).head
    assert(homGt.getVariant.getStart === 839355L)
    assert(homGt.getVariant.getEnd === 839356L)
    assert(homGt.getVariant.getReferenceAllele === "A")
    assert(homGt.getVariant.getAlternateAllele === "C")
  }
}
