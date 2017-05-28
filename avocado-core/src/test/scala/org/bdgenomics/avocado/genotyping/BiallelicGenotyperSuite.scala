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

import org.bdgenomics.adam.models.{
  RecordGroup,
  RecordGroupDictionary,
  SequenceDictionary,
  SequenceRecord
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
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
  lazy val summaryObservations = ScoredObservation.createScores(sc, 90, 90, 2)
    .collect()
    .toSeq

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

  test("properly handle haploid genotype state") {
    val (state, qual) = BiallelicGenotyper.genotypeStateAndQuality(
      Array(5.0, 0.394829))

    assert(state === 0)
    assert(MathUtils.fpEquals(qual, 20.0, tol = 1e-3))
  }

  test("properly handle diploid genotype state with het call") {
    val (state, qual) = BiallelicGenotyper.genotypeStateAndQuality(
      Array(-10.0, 5.0, -1.907755))

    assert(state === 1)
    assert(MathUtils.fpEquals(qual, 30.0, tol = 1e-3))
  }

  test("properly handle triploid genotype state with hom alt call") {
    val (state, qual) = BiallelicGenotyper.genotypeStateAndQuality(
      Array(-10.0, 5.0, -1.907755, 14.210340))

    assert(state === 3)
    assert(MathUtils.fpEquals(qual, 40.0, tol = 1e-3))
  }

  test("scoring read that overlaps no variants should return empty observations") {
    val perfectScores = BiallelicGenotyper.readToObservations((perfectRead, Iterable.empty), 2)
    assert(perfectScores.isEmpty)

    val snpScores = BiallelicGenotyper.readToObservations((snpRead, Iterable.empty), 2)
    assert(snpScores.isEmpty)
  }

  sparkTest("score snp in a read with no evidence of the snp") {
    val scores = BiallelicGenotyper.readToObservations(
      (perfectRead, Iterable(snp)), 2)
    assert(scores.size === 1)

    val (snpVariant, snpSumObservation) = scores.head
    val snpObservation = snpSumObservation.toObservation(summaryObservations)
    assert(snpVariant.toVariant === snp)
    assert(snpObservation.squareMapQ === 50 * 50)
    assert(snpObservation.alleleCoverage === 0)
    assert(snpObservation.otherCoverage === 1)
    assert(snpObservation.alleleForwardStrand === 0)
    assert(snpObservation.otherForwardStrand === 1)
    assert(snpObservation.copyNumber === 2)
    assert(MathUtils.fpEquals(snpObservation.referenceLogLikelihoods(0), os.logL(0, 2, 0.9999, 0.99999)))
    assert(MathUtils.fpEquals(snpObservation.referenceLogLikelihoods(1), os.logL(1, 2, 0.9999, 0.99999)))
    assert(MathUtils.fpEquals(snpObservation.referenceLogLikelihoods(2), os.logL(2, 2, 0.9999, 0.99999)))
  }

  sparkTest("score snp in a read with evidence of the snp") {
    val scores = BiallelicGenotyper.readToObservations(
      (snpRead, Iterable(snp)), 2)
    assert(scores.size === 1)

    val (snpVariant, snpSumObservation) = scores.head
    val snpObservation = snpSumObservation.toObservation(summaryObservations)
    assert(snpVariant.toVariant === snp)
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
      Array(-24.0, -4.0, -12.0),
      Array(-10.0, -1.0, -10.0),
      Array(0.0, 0.0, 0.0),
      7, 9)
    val genotype = BiallelicGenotyper.observationToGenotype((snp, obs), "sample")

    assert(genotype.getVariant === snp)
    assert(genotype.getStart === snp.getStart)
    assert(genotype.getEnd === snp.getEnd)
    assert(genotype.getContigName === snp.getContigName)
    assert(genotype.getSampleId === "sample")
    assert(genotype.getGenotypeQuality === 73)
    assert(genotype.getAlleles.size === 2)
    assert(genotype.getAlleles.get(0) === GenotypeAllele.ALT)
    assert(genotype.getAlleles.get(1) === GenotypeAllele.REF)
    assert(genotype.getReferenceReadDepth === 9)
    assert(genotype.getAlternateReadDepth === 7)
    assert(genotype.getStrandBiasComponents.size === 4)
    assert(genotype.getStrandBiasComponents.get(0) === 5)
    assert(genotype.getStrandBiasComponents.get(1) === 4)
    assert(genotype.getStrandBiasComponents.get(2) === 3)
    assert(genotype.getStrandBiasComponents.get(3) === 4)
    assert(genotype.getGenotypeLikelihoods.size === 3)
  }

  ignore("force call possible STR/indel") {
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
      !gt.getAlleles.forall(_ == GenotypeAllele.REF)
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
            !gt.getAlleles.forall(_ == GenotypeAllele.REF) &&
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
    assert(gt.getAlleles.count(_ == GenotypeAllele.REF) === 1)
    assert(gt.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
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
      TestHardFilterGenotypesArgs(minHetIndelAltAllelicFraction = -1f,
        maxHetIndelAltAllelicFraction = -1f,
        minHomIndelAltAllelicFraction = -1f))
      .transform(rdd => {
        rdd.filter(gt => {
          gt.getVariantCallingAnnotations
            .getFiltersPassed && gt.getAlternateReadDepth != 0
        })
      })

    val gts = filteredGenotypes.rdd.collect

    assert(gts.size === 1)
    val gt = gts.head
    assert(gt.getVariant.getStart === 832735L)
    assert(gt.getVariant.getEnd === 832736L)
    assert(gt.getVariant.getReferenceAllele === "A")
    assert(gt.getVariant.getAlternateAllele === "AGTTTT")
    assert(gt.getAlleles.count(_ == GenotypeAllele.REF) === 1)
    assert(gt.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
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
            .getFiltersPassed &&
            gt.getReadDepth > 30
        })
      })

    val gts = filteredGenotypes.rdd.collect
    //assert(gts.size === 2)

    /*
    val hetGt = gts.filter(gt => {
      gt.getAlleles.count(_ == GenotypeAllele.Alt) == 1
    }).head
    assert(hetGt.getVariant.getStart === 839404L)
    assert(hetGt.getVariant.getEnd === 839405L)
    assert(hetGt.getVariant.getReferenceAllele === "G")
    assert(hetGt.getVariant.getAlternateAllele === "A")
    */
    val homGt = gts.filter(gt => {
      gt.getAlleles.count(_ == GenotypeAllele.ALT) == 2
    }).head
    assert(homGt.getVariant.getStart === 839355L)
    assert(homGt.getVariant.getEnd === 839356L)
    assert(homGt.getVariant.getReferenceAllele === "A")
    assert(homGt.getVariant.getAlternateAllele === "C")
  }

  sparkTest("score a single read covering a deletion") {
    val readPath = resourceUrl("NA12878.chr1.567239.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getReadName == "H06JUADXX130110:1:2102:2756:44620")
          .cache
      })
    assert(reads.rdd.count === 1)

    val variants = DiscoverVariants(reads)
      .transform(rdd => {
        rdd.filter(v => {
          v.getStart == 567238L &&
            v.getEnd == 567240L &&
            v.getReferenceAllele == "CG" &&
            v.getAlternateAllele == "C"
        })
      }).rdd.collect
    assert(variants.size === 1)

    val obs = BiallelicGenotyper.readToObservations((reads.rdd.first,
      variants.toIterable), 2)

    assert(obs.size === 1)
  }

  sparkTest("discover and force call hom alt deletion") {
    val readPath = resourceUrl("NA12878.chr1.567239.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val variants = DiscoverVariants(reads)
      .transform(rdd => {
        rdd.filter(v => {
          v.getStart == 567238L &&
            v.getEnd == 567240L &&
            v.getReferenceAllele == "CG" &&
            v.getAlternateAllele == "C"
        }).cache
      })

    assert(variants.rdd.count === 1)

    val gts = BiallelicGenotyper.call(reads, variants, 2)
    val gtArray = gts.rdd.collect
    assert(gtArray.size === 1)
    val gt = gtArray.head

    assert(gt.getAlleles.forall(_ == GenotypeAllele.ALT))
  }

  sparkTest("call hom alt AGCCAGTGGACGCCGACCT->A deletion at 1/875159") {
    val readPath = resourceUrl("NA12878.chr1.875159.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val variants = DiscoverVariants(reads)
      .transform(rdd => {
        rdd.filter(v => {
          v.getStart == 875158L &&
            v.getEnd == 875177L &&
            v.getReferenceAllele == "AGCCAGTGGACGCCGACCT" &&
            v.getAlternateAllele == "A"
        }).cache
      })
    assert(variants.rdd.count === 1)

    val gts = BiallelicGenotyper.call(reads, variants, 2)
    val gtArray = gts.rdd.collect
    assert(gtArray.size === 1)
    val gt = gtArray.head

    assert(gt.getAlleles.forall(_ == GenotypeAllele.ALT))
  }

  sparkTest("call hom alt C->G snp at 1/877715") {
    val readPath = resourceUrl("NA12878.chr1.877715.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val gts = BiallelicGenotyper.discoverAndCall(reads, 2)
      .transform(rdd => {
        rdd.filter(gt => gt.getStart == 877714)
      })
    val gtArray = gts.rdd.collect
    assert(gtArray.size === 1)
    val gt = gtArray.head

    assert(gt.getAlleles.forall(_ == GenotypeAllele.ALT))
  }

  sparkTest("call hom alt ACAG->A deletion at 1/886049") {
    val readPath = resourceUrl("NA12878.chr1.886049.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val gts = BiallelicGenotyper.discoverAndCall(reads, 2)
      .transform(rdd => {
        rdd.filter(gt => gt.getStart == 886048)
      })
    val gtArray = gts.rdd.collect
    assert(gtArray.size === 1)
    val gt = gtArray.head

    assert(gt.getAlleles.forall(_ == GenotypeAllele.ALT))
  }

  sparkTest("call hom alt GA->CC mnp at 1/889158–9") {
    val readPath = resourceUrl("NA12878.chr1.889159.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val gts = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optPhredThreshold = Some(30)).transform(rdd => {
        rdd.filter(gt => !gt.getAlleles.forall(_ == GenotypeAllele.REF))
      })
    val gtArray = gts.rdd.collect
    assert(gtArray.size === 2)
    assert(gtArray.count(gt => {
      gt.getVariant.getStart == 889157L &&
        gt.getVariant.getEnd == 889158L &&
        gt.getVariant.getReferenceAllele == "G" &&
        gt.getVariant.getAlternateAllele == "C"
    }) === 1)
    assert(gtArray.count(gt => {
      gt.getVariant.getStart == 889158L &&
        gt.getVariant.getEnd == 889159L &&
        gt.getVariant.getReferenceAllele == "A" &&
        gt.getVariant.getAlternateAllele == "C"
    }) === 1)
    gtArray.foreach(gt => {
      assert(gt.getAlleles.forall(_ == GenotypeAllele.ALT))
    })
  }

  sparkTest("call hom alt C->CCCCT insertion at 1/866511") {
    val readPath = resourceUrl("NA12878.chr1.866511.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      }).realignIndels()

    val gts = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optPhredThreshold = Some(30)).transform(rdd => {
        rdd.filter(gt => gt.getAlleles.forall(_ == GenotypeAllele.ALT))
      })
  }

  sparkTest("call het ATG->A deletion at 1/905130") {
    val readPath = resourceUrl("NA12878.chr1.905130.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val gts = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optPhredThreshold = Some(30)).transform(rdd => {
        rdd.filter(gt => {
          !gt.getAlleles.forall(_ == GenotypeAllele.REF) &&
            gt.getVariant.getStart == 905129L
        })
      }).rdd.collect

    assert(gts.size === 1)
    val gt = gts.head
    assert(gt.getVariant.getContigName === "1")
    assert(gt.getVariant.getStart === 905129L)
    assert(gt.getVariant.getEnd === 905132L)
    assert(gt.getVariant.getReferenceAllele === "ATG")
    assert(gt.getVariant.getAlternateAllele === "A")
    assert(gt.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
  }

  sparkTest("call het AG->A deletion at 1/907170") {
    val readPath = resourceUrl("NA12878.chr1.907170.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 0)
      })

    val gts = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optPhredThreshold = Some(30)).transform(rdd => {
        rdd.filter(gt => !gt.getAlleles.forall(_ == GenotypeAllele.REF))
      }).rdd.collect

    assert(gts.size === 1)
    val gt = gts.head
    assert(gt.getVariant.getContigName === "1")
    assert(gt.getVariant.getStart === 907169L)
    assert(gt.getVariant.getEnd === 907171L)
    assert(gt.getVariant.getReferenceAllele === "AG")
    assert(gt.getVariant.getAlternateAllele === "A")
    assert(gt.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
  }

  sparkTest("call het T->G snp at 1/240898") {
    val readPath = resourceUrl("NA12878.chr1.240898.sam")
    val reads = sc.loadAlignments(readPath.toString)
      .transform(rdd => {
        rdd.filter(_.getMapq > 10)
      })

    val gts = BiallelicGenotyper.discoverAndCall(reads,
      2,
      optPhredThreshold = Some(25)).transform(rdd => {
        rdd.filter(gt => !gt.getAlleles.forall(_ == GenotypeAllele.REF))
      }).rdd.collect

    assert(gts.size === 1)
    val gt = gts.head
    assert(gt.getVariant.getContigName === "1")
    assert(gt.getVariant.getStart === 240897L)
    assert(gt.getVariant.getEnd === 240898L)
    assert(gt.getVariant.getReferenceAllele === "T")
    assert(gt.getVariant.getAlternateAllele === "G")
    assert(gt.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
  }

  sparkTest("make het alt calls at biallelic snp locus") {
    def makeRead(allele: Char): AlignmentRecord = {
      assert(allele != 'T')
      AlignmentRecord.newBuilder
        .setContigName("ctg")
        .setStart(10L)
        .setEnd(15L)
        .setSequence("AC%sTG".format(allele))
        .setCigar("5M")
        .setMismatchingPositions("2T2")
        .setQual(Seq(50, 50, 50, 50, 50).map(q => (q + 33).toInt).mkString)
        .setMapq(50)
        .setReadMapped(true)
        .setPrimaryAlignment(true)
        .build
    }

    val reads = Seq(makeRead('A'), makeRead('A'), makeRead('A'), makeRead('A'),
      makeRead('C'), makeRead('C'), makeRead('C'), makeRead('C'))
    val readRdd = AlignmentRecordRDD(
      sc.parallelize(reads),
      SequenceDictionary(SequenceRecord("ctg", 16L)),
      RecordGroupDictionary(Seq(RecordGroup("rg1", "rg1"))))

    val gts = BiallelicGenotyper.discoverAndCall(readRdd,
      2).rdd.collect
    assert(gts.size === 2)
    assert(gts.forall(gt => gt.getAlleles.size == 2))
    assert(gts.forall(gt => gt.getAlleles.count(_ == GenotypeAllele.ALT) == 1))
    assert(gts.forall(gt => gt.getAlleles.count(_ == GenotypeAllele.OTHER_ALT) == 1))
    val alleles = gts.map(gt => gt.getVariant.getAlternateAllele).toSet
    assert(alleles("A"))
    assert(alleles("C"))
  }
}
