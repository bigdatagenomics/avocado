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

package edu.berkeley.cs.amplab.avocado.calls.pileup

import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup, Base, ADAMGenotype, VariantType}
import edu.berkeley.cs.amplab.adam.models.{ADAMRod, ADAMVariantContext}
import edu.berkeley.cs.amplab.adam.util.PhredUtils
import edu.berkeley.cs.amplab.avocado.calls.VariantCallCompanion
import edu.berkeley.cs.amplab.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import scala.math.pow
import scala.collection.JavaConversions._

object PileupCallSimpleSNP extends VariantCallCompanion {

  val callName = "SimpleSNP"

  def apply (stats: AvocadoConfigAndStats,
             config: SubnodeConfiguration): PileupCallSimpleSNP = {

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
class PileupCallSimpleSNP (ploidy: Int) extends PileupCall {
  
  val companion: VariantCallCompanion = PileupCallSimpleSNP

  /**
   * Takes pileup info and likelihoods and writes out to a variant list.
   *
   * @param pileupHead Single pileup that contains location and sample info.
   * @param likelihood List of likelihoods for homozygous ref/heterozygous/homozygous non ref.
   * @param maxNonRefBase Highest likelihood base for mismatch.
   * @return List of variants called, with 0 (homozygous reference) or 1 entry.
   */
  protected def writeCallInfo (pileupHead: ADAMPileup, 
                               likelihood: List[Double], 
                               maxNonRefBase: Option[Base]): List[ADAMGenotype] = {
    assert(likelihood.length == 3)

    // get phred scores
    val homozygousRefPhred = PhredUtils.successProbabilityToPhred(likelihood(0))
    val heterozygousPhred = PhredUtils.successProbabilityToPhred(likelihood(1))
    val homozygousNonPhred = PhredUtils.successProbabilityToPhred(likelihood(2))

    // simplifying assumption - snps are biallelic
    val call = if (likelihood.indexOf(likelihood.max) != 0 &&
                   maxNonRefBase.isEmpty) {
      // genotype likelihood is not in favor of reference, but we do not have any non ref bases?
      log.warn("Want to call non-homozygous ref genotype, but don't see non-ref bases @ " +
               pileupHead.getPosition + ".")
      List[ADAMGenotype]()      
    } else if (likelihood.indexOf(likelihood.max) == 1) {
      // hetereozygous
      assert(pileupHead.getReferenceBase != maxNonRefBase.get,
             "Cannot have reference be equal to non-ref base, @ " + pileupHead.getPosition + ".")
      
      val genotypeRef = ADAMGenotype.newBuilder()
        .setReferenceId(pileupHead.getReferenceId)
        .setReferenceName(pileupHead.getReferenceName)
        .setPosition(pileupHead.getPosition)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setAlleleVariantType(VariantType.SNP)
        .setGenotypeQuality(heterozygousPhred)
        .setPloidy(2)
        .setHaplotypeNumber(0)
        .setIsReference(true)
        .setReferenceAllele(pileupHead.getReferenceBase.toString)
        .setAllele(pileupHead.getReferenceBase.toString)
        .setExpectedAlleleDosage(1.0)
	.build()
      val genotypeNonRef = ADAMGenotype.newBuilder()
        .setReferenceId(pileupHead.getReferenceId)
        .setReferenceName(pileupHead.getReferenceName)
        .setPosition(pileupHead.getPosition)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setAlleleVariantType(VariantType.SNP)
        .setGenotypeQuality(heterozygousPhred)
        .setPloidy(2)
        .setHaplotypeNumber(1)
        .setIsReference(false)
        .setReferenceAllele(pileupHead.getReferenceBase.toString)
        .setAllele(maxNonRefBase.get.toString)
        .setExpectedAlleleDosage(1.0)
	.build()

      List(genotypeRef, genotypeNonRef)
    } else if (likelihood.indexOf(likelihood.max) == 2) {
      // homozygous non reference
      assert(pileupHead.getReferenceBase != maxNonRefBase.get,
             "Cannot have reference be equal to non-ref base, @ " + pileupHead.getPosition + ".")
            
      val genotypeNonRef0 = ADAMGenotype.newBuilder()
        .setReferenceId(pileupHead.getReferenceId)
        .setReferenceName(pileupHead.getReferenceName)
        .setPosition(pileupHead.getPosition)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setAlleleVariantType(VariantType.SNP)
        .setGenotypeQuality(homozygousNonPhred)
        .setPloidy(2)
        .setHaplotypeNumber(0)
        .setIsReference(false)
        .setReferenceAllele(pileupHead.getReferenceBase.toString)
        .setAllele(maxNonRefBase.get.toString)
        .setExpectedAlleleDosage(2.0)
	.build()
      val genotypeNonRef1 = ADAMGenotype.newBuilder()
        .setReferenceId(pileupHead.getReferenceId)
        .setReferenceName(pileupHead.getReferenceName)
        .setPosition(pileupHead.getPosition)
        .setSampleId(pileupHead.getRecordGroupSample)
        .setAlleleVariantType(VariantType.SNP)
        .setGenotypeQuality(homozygousNonPhred)
        .setPloidy(2)
        .setHaplotypeNumber(1)
        .setIsReference(false)
        .setReferenceAllele(pileupHead.getReferenceBase.toString)
        .setAllele(maxNonRefBase.get.toString)
        .setExpectedAlleleDosage(2.0)
	.build()

      List(genotypeNonRef0, genotypeNonRef1)
    } else {
      // homozygous
      List[ADAMGenotype]()
    }

    if (call.length != 0) {
      assert(call.length == ploidy,
             "Calls emitted must be consistent with ploidy.")
    }

    call
  }

  /**
   * Scores likelihoods for genotypes using equation from Li 2009.
   *
   * @param[in] pileup Rod containing pileup bases.
   * @return List of doubles corresponding to likelihood of homozygous (0), heterozygous (1), and homozygous non-reference (2).
   */
  def scoreGenotypeLikelihoods (pileup: List[ADAMPileup]): List[Double] = {
    
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
     * - g in 0..2 is genotype --> homozygous reference, hetero, other
     * - e is error probablity of nucleotide
     * - m is ploidy
     * - k is number of bases that match reference
     * - l is the number of bases that mismatch against reference
     */
    val likelihoods = states.map(g => {
      // contribution of bases that match the reference
      val refBases = pileup.filter(v => v.getReadBase == v.getReferenceBase)
      val productMatch = if (refBases.length != 0) {
        refBases.map((base) => {
          val epsilon = 1.0 - PhredUtils.phredToSuccessProbability((base.getMapQuality + base.getSangerQuality) / 2)
        
	  pow((ploidy - g) * epsilon + g * (1 - epsilon), base.getCountAtPosition.toDouble)
        }).reduce(_ * _)
      } else {
        1.0
      }

      // contribution of bases that do not match the base
      val mismatchBases = pileup.filter(v => {
	val readBase = Option(v.getReadBase) match {
	  case Some(base) => base
	  case None => Base.N // TODO: add better exception handling code - this case shouldn't happen unless there is deletion evidence in the rod
	}

	((readBase != v.getReferenceBase) && v.getRangeLength == 0)
      })

      val productMismatch = if (mismatchBases.length != 0) {
        mismatchBases.map((base) => {
          val epsilon = 1.0 - PhredUtils.phredToSuccessProbability((base.getMapQuality + base.getSangerQuality) / 2)
          
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
    def pickMaxBase (kv1: (Base, Int), kv2: (Base, Int)): (Base, Int) = {
      if (kv1._2 > kv2._2) {
	kv1
      } else {
	kv2
      }
    }
    
    // reduce down to get the base with the highest count
    if (nonRefBaseCount.isEmpty){
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
   * @param[in] pileup List of pileups. Should only contain one rod.
   * @return List of variants seen at site. List can contain 0 or 1 elements - value goes to flatMap.
   */
  protected def callSNP (pileup: List[ADAMPileup]): List[ADAMVariantContext] = {
    
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
   * @param[in] pileupGroups An RDD containing lists of pileups.
   * @return An RDD containing called variants.
   */
  override def callRods (pileups: RDD [ADAMRod]): RDD [ADAMVariantContext] = {
    log.info("Calling SNPs on pileups and flattening.")
    pileups.map(_.pileups)
      .map(callSNP)
      .flatMap((p: List[ADAMVariantContext]) => p)
  }

  override def isCallable (): Boolean = true
}


