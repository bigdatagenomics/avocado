/*
 * Copyright (c) 2013. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.avocado.calls.pileup

import org.bdgenomics.adam.avro.{ Base, ADAMContig, ADAMGenotype, ADAMPileup, ADAMVariant }
import org.bdgenomics.adam.models.{ ADAMRod, ADAMVariantContext }
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.calls.VariantCallCompanion
import org.bdgenomics.avocado.partitioners.PartitionSet
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import scala.math.pow
import scala.collection.JavaConversions._

object PileupCallSimpleSNP extends VariantCallCompanion {

  val callName = "SimpleSNP"

  def apply(stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration,
            partitions: PartitionSet): PileupCallSimpleSNP = {

    val ploidy = config.getInt("ploidy", 2)

    new PileupCallSimpleSNP(ploidy)
  }
}

/**
 * Class to call SNPs. Implements the SNP genotyping methods described in:
 *
 * Li, Heng. "A statistical framework for SNP calling, mutation discovery, association mapping and population
 * genetical parameter estimation from sequencing data." Bioinformatics 27.21 (2011): 2987-2993.
 *
 * for a single sample. As we are taking in a single sample, we do not calculate a minor allele frequency (MAF).
 * We only call the genotype if the likelihood has a phred score greater than or equal to 30.
 * At the current point in time, we assume that we are running on a diploid organism.
 */
class PileupCallSimpleSNP(ploidy: Int) extends PileupCall {

  val companion: VariantCallCompanion = PileupCallSimpleSNP

  /**
   * Takes pileup info and likelihoods and writes out to a variant list.
   * TODO: under the new schema we don't need to return a list of
   * ADAMGenotype, since each ADAMGenotype contains the same
   * ADAMVariant.  Figure out what the best return interface is.
   *
   * @param pileupHead Single pileup that contains location and sample info.
   * @param likelihood List of likelihoods for homozygous ref/heterozygous/homozygous non ref.
   * @param maxNonRefBase Highest likelihood base for mismatch.
   * @return List of variants called, with 0 (homozygous reference) or 1 entry.
   */
  protected def writeCallInfo(pileupHead: ADAMPileup,
                              likelihood: List[Double],
                              maxNonRefBase: Option[Base]): List[ADAMGenotype] = {
    assert(likelihood.length == 3)

    // get phred scores
    val homozygousRefPhred = PhredUtils.successProbabilityToPhred(likelihood(2))
    val heterozygousPhred = PhredUtils.successProbabilityToPhred(likelihood(1))
    val homozygousNonPhred = PhredUtils.successProbabilityToPhred(likelihood(0))

    // simplifying assumption - snps are biallelic
    val call = if (likelihood.indexOf(likelihood.max) != 2 &&
      maxNonRefBase.isEmpty) {
      // genotype likelihood is not in favor of reference, but we do not have any non ref bases?
      log.warn("Want to call non-homozygous ref genotype, but don't see non-ref bases @ " +
        pileupHead.getPosition + ".")
      List[ADAMGenotype]()
    } else if (likelihood.indexOf(likelihood.max) == 1) {
      // hetereozygous
      assert(pileupHead.getReferenceBase != maxNonRefBase.get,
        "Cannot have reference be equal to non-ref base, @ " + pileupHead.getPosition + ".")

      val contig = ADAMContig.newBuilder
        .setContigName(pileupHead.getContig.getContigName)
        .build
      val variant = ADAMVariant.newBuilder
        .setContig(contig)
        .setPosition(pileupHead.getPosition)
        .setReferenceAllele(pileupHead.getReferenceBase.toString)
        .setVariantAllele(maxNonRefBase.get.toString)
        .build
      val genotypeRef = ADAMGenotype.newBuilder()
        .setVariant(variant)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setGenotypeQuality(heterozygousPhred)
        .setExpectedAlleleDosage(1.0f)
        .build()
      val genotypeNonRef = ADAMGenotype.newBuilder()
        .setVariant(variant)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setGenotypeQuality(heterozygousPhred)
        .setExpectedAlleleDosage(1.0f)
        .build()

      List(genotypeRef, genotypeNonRef)
    } else if (likelihood.indexOf(likelihood.max) == 0) {
      // homozygous non reference
      assert(pileupHead.getReferenceBase != maxNonRefBase.get,
        "Cannot have reference be equal to non-ref base, @ " + pileupHead.getPosition + ".")

      val contig = ADAMContig.newBuilder
        .setContigName(pileupHead.getContig.getContigName)
        .build
      val variant = ADAMVariant.newBuilder
        .setContig(contig)
        .setPosition(pileupHead.getPosition)
        .setReferenceAllele(pileupHead.getReferenceBase.toString)
        .setVariantAllele(maxNonRefBase.get.toString)
        .build
      val genotypeNonRef0 = ADAMGenotype.newBuilder()
        .setVariant(variant)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setGenotypeQuality(homozygousNonPhred)
        .setExpectedAlleleDosage(2.0f)
        .build()
      val genotypeNonRef1 = ADAMGenotype.newBuilder()
        .setVariant(variant)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setGenotypeQuality(homozygousNonPhred)
        .setExpectedAlleleDosage(2.0f)
        .build()

      List(genotypeNonRef0, genotypeNonRef1)
    } else {
      // homozygous
      List[ADAMGenotype]()
    }

    if (call.length != 0) {
      log.info("Writing call info at " + pileupHead.getPosition)
    }

    call
  }

  /**
   * Scores likelihoods for genotypes using equation from Li 2009.
   *
   * @param[in] pileup Rod containing pileup bases.
   * @return List of doubles corresponding to likelihood of homozygous (2), heterozygous (1), and homozygous non-reference (0).
   */
  def scoreGenotypeLikelihoods(pileup: List[ADAMPileup]): List[Double] = {

    // count bases in pileup
    val k = pileup.map(_.getCountAtPosition).reduce(_ + _)

    // genotype states
    val states = List(0, 1, 2)

    /* find genotype for pileup
     * likelihood for genotype is derived from:
     * L(g) = 1/m^k *
     *        product over j->1..l ((m - g) * e + g * (1 - e)) *
     *        product over j->l+1..k ((m - g) * (1 - e) + g * e)
     *
     * Where:
     * - g in 0..2 is genotype -> g ==  # of reference bases in allele
     * - e is error probablity of nucleotide
     * - m is ploidy
     * - k is number of bases that match reference
     * - l is the number of bases that mismatch against reference
     */
    val refBases = pileup.filter(v => v.getReadBase == v.getReferenceBase)
    val mismatchBases = pileup.filter(v => {
      val readBase = Option(v.getReadBase) match {
        case Some(base) => base
        case None       => Base.N // TODO: add better exception handling code - this case shouldn't happen unless there is deletion evidence in the rod
      }

      ((readBase != v.getReferenceBase) && (v.getRangeLength == 0 || v.getRangeLength == null))
    })

    val likelihoods = states.map(g => {
      // contribution of bases that match the reference
      val productMatch = if (refBases.length != 0) {
        refBases.map((base) => {
          val epsilon = PhredUtils.phredToErrorProbability((base.getMapQuality + base.getSangerQuality) / 2)

          pow((ploidy - g) * epsilon + g * (1 - epsilon), base.getCountAtPosition.toDouble)
        }).reduce(_ * _)
      } else {
        1.0
      }

      // contribution of bases that do not match the base
      val productMismatch = if (mismatchBases.length != 0) {
        mismatchBases.map((base) => {
          val epsilon = PhredUtils.phredToErrorProbability((base.getMapQuality + base.getSangerQuality) / 2)

          pow((ploidy - g) * (1 - epsilon) + g * epsilon, base.getCountAtPosition.toDouble)
        }).reduce(_ * _)
      } else {
        1.0
      }

      productMatch * productMismatch / pow(ploidy.toDouble, k.toDouble)
    })

    likelihoods
  }

  /**
   * Picks the most frequently seen non-reference base. If no non-reference bases are
   * seen, then emits None.
   *
   * @param pileup List of pileups. Assume that all pileups are at a single locus, and do not
   * show SNP evidence.
   * @return Option containing most recently seen non-reference base if found, else none.
   */
  def getMaxNonRefBase(pileup: List[ADAMPileup]): Option[Base] = {

    // get a count of the total number of bases that are a mismatch with the reference
    val nonRefBaseCount = pileup.filter(r => r.getReadBase != r.getReferenceBase)
      .groupBy(_.getReadBase)
      .map(kv => (kv._1, kv._2.length))

    /**
     * Out of two Base, Int pairs, picks the one with the higher count.
     *
     * @param[in] kv1 Key/value pair containing a Base and it's count.
     * @param[in] kv2 Key/value pair containing a Base and it's count.
     * @return The key/value pair with the higher count.
     */
    def pickMaxBase(kv1: (Base, Int), kv2: (Base, Int)): (Base, Int) = {
      if (kv1._2 > kv2._2) {
        kv1
      } else {
        kv2
      }
    }

    // reduce down to get the base with the highest count
    if (nonRefBaseCount.isEmpty) {
      None
    } else {
      Some(nonRefBaseCount.reduce(pickMaxBase)._1)
    }
  }

  /**
   * Compensates likelihoods.
   *
   * @param likelihoods List of likelihood values.
   * @param pileup List of pileups at this position.
   * @return Likelihoods after compensation.
   */
  def compensate(likelihoods: List[Double], pileup: List[ADAMPileup]): List[Double] = {
    likelihoods
  }

  /**
   * Calls a SNP for a single pileup, if there is sufficient evidence of a mismatch from the reference.
   * Takes in a single pileup rod from a single sample at a single locus. For simplicity, we assume that
   * all sites are biallelic.
   *
   * @param[in] rod ADAMRod
   * @return List of variants seen at site. List can contain 0 or 1 elements - value goes to flatMap.
   */
  protected def callSNP(rod: ADAMRod): List[ADAMVariantContext] = {
    val pileup = rod.pileups
    if (pileup.forall(_.getRangeLength == null)) {
      val loci = pileup.head.getPosition
      log.info("Calling pileup at " + loci)

      // score off of rod info
      val likelihood = scoreGenotypeLikelihoods(pileup)

      // compensate
      val compensatedLikelihood = compensate(likelihood, pileup)

      // get most often seen non-reference base
      val maxNonRefBase = getMaxNonRefBase(pileup)

      // write calls to list
      genotypesToVariantContext(writeCallInfo(pileup.head, likelihood, maxNonRefBase))
    } else {
      // cannot call indels
      log.warn("Saw indel evidence at " + pileup.head.getPosition + ". Ignored.")
      List[ADAMVariantContext]()
    }
  }

  /**
   * Call variants using simple pileup based SNP calling algorithm.
   *
   * @param[in] pileups An RDD containing of ADAMRods.
   * @return An RDD containing called variants.
   */
  override def callRods(pileups: RDD[ADAMRod]): RDD[ADAMVariantContext] = {
    log.info("Calling SNPs on pileups and flattening.")
    pileups
      .map(callSNP)
      .flatMap((p: List[ADAMVariantContext]) => p)
  }

  override def isCallable(): Boolean = true
}

