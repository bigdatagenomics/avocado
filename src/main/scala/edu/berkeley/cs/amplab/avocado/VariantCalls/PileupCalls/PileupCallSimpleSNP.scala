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

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup, ADAMVariant, Base, ADAMGenotype, VariantType}
import edu.berkeley.cs.amplab.adam.models.ADAMRod
import edu.berkeley.cs.amplab.avocado.utils.Phred
import scala.math.pow
import scala.collection.mutable.MutableList
import scala.collection.JavaConversions._

/**
 * Class to call SNPs. Implements the SNP genotyping methods described in:
 *
 * Li, Heng. "A statistical framework for SNP calling, mutation discovery, association mapping and population
 * genetical parameter estimation from sequencing data." Bioinformatics 27.21 (2011): 2987-2993.
 *
 * for a single sample. As we are taking in a single sample, we do not calculate a minor allele frequency (MAF); rather,
 * we assume a MAF of 1 in 10. We only call the genotype if the likelihood has a phred score greater than or equal to 30.
 * At the current point in time, we assume that we are running on a diploid organism.
 */
class PileupCallSimpleSNP extends PileupCall {
  
  val callName = "SimpleSNP"

  // at current point in time, assume human diploid
  // TODO: extend to arbitrary ploidy
  val ploidy = 2

  // minor allele frequency - currently set to 1 in 10
  val minorAlleleFrequency = 0.1

  // PHRED score threshold for calling
  val variantQuality = 30

  /**
   * Takes pileup info and likelihoods and writes out to a variant list.
   *
   * @param[in] pileupHead Single pileup that contains location and sample info.
   * @param[in] likelihood List of likelihoods for homozygous ref/heterozygous/homozygous non ref.
   * @param[in] maxNonRefBase Highest likelihood base for mismatch.
   * @return List of variants called, with 0 (homozygous reference) or 1 entry.
   */
  protected def writeCallInfo (pileupHead: ADAMPileup, likelihood: MutableList[Double], maxNonRefBase: Base): (List[ADAMVariant], List[ADAMGenotype]) = {
    assert (likelihood.length == 3)

    log.info ("Writing call info at " + pileupHead.getPosition)

    // get phred scores
    val homozygousRefPhred = Phred.probabilityToPhred (likelihood (0))
    val heterozygousPhred = Phred.probabilityToPhred (likelihood (1))
    val homozygousNonPhred = Phred.probabilityToPhred (likelihood (2))

    // build genotype and variant descriptors and write out
    // TODO: some code here converts List[Int] to java.lang.List[java.lang.Integer] and is ugly
    // add implicit conversion to utils

    // simplifying assumption - snps are biallelic
    if (likelihood.indexOf (likelihood.max) == 1 &&
	heterozygousPhred >= variantQuality) {
      // hetereozygous
      // and, variant quality is >= threshold 
            
      val genotype = ADAMGenotype.newBuilder ()
	.setSampleId (pileupHead.getRecordGroupSample)
	.setGenotype ("0,1")
	.setPhredLikelihoods (heterozygousPhred.toString)
	.build()

      val variant = ADAMVariant.newBuilder ()
	.setReferenceName (pileupHead.getReferenceName)
	.setStartPosition (pileupHead.getPosition)
	.setReferenceAllele (pileupHead.getReferenceBase.toString)
	.setAlternateAlleles (maxNonRefBase.toString)
	.setAlleleCount (2)
	.setType (VariantType.SNP)
	.build()

      (List (variant), List (genotype))
    } else if (likelihood.indexOf (likelihood.max) == 2 &&
	       homozygousNonPhred >= variantQuality) {
      // homozygous non reference
      // and, variant quality is >= threshold
      
      val genotype = ADAMGenotype.newBuilder ()
	.setSampleId (pileupHead.getRecordGroupSample)
	.setGenotype ("1,1")
	.setPhredLikelihoods (homozygousNonPhred.toString)
	.build()

      val variant = ADAMVariant.newBuilder ()
	.setReferenceName (pileupHead.getReferenceName)
	.setStartPosition (pileupHead.getPosition)
	.setReferenceAllele (pileupHead.getReferenceBase.toString)
	.setAlternateAlleles (maxNonRefBase.toString)
	.setAlleleCount (1)
	.setType (VariantType.SNP)
	.build()

      (List (variant), List (genotype))
    } else {
      // homozygous
      // or, variant quality is < threshold
      (List [ADAMVariant] (), List [ADAMGenotype] ())
    }
  }

  /**
   * Scores likelihoods for genotypes using equation from Li 2009.
   *
   * @param[in] pileup Rod containing pileup bases.
   * @return List of doubles corresponding to likelihood of homozygous (0), heterozygous (1), and homozygous non-reference (2).
   */
  protected def scoreGenotypeLikelihoods (pileup: List[ADAMPileup]): MutableList[Double] = {

    // allocate list for homozygous reference, heterozygous, homozygous minor allele
    var likelihood = MutableList[Double](0.0, 0.0, 0.0)

    // count bases in pileup
    val k = pileup.map (_.getCountAtPosition).reduce (_ + _)

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
    for (g <- 0 to 2) {
      // contribution of bases that match the reference
      val refBases = pileup.filter(v => v.getReadBase == v.getReferenceBase)
      val productMatch = if (refBases.length != 0) {
        refBases.map ((base) => {
          val epsilon = Phred.phredToProbability ((base.getMapQuality + base.getSangerQuality) / 2)
        
	  pow ((ploidy - g) * epsilon + g * (1 - epsilon), base.getCountAtPosition.toDouble)
        }).reduce (_ * _)
      } else {
        0.0
      }

      // contribution of bases that do not match the base
      val mismatchBases = pileup.filter(v => {
	val readBase = Option (v.getReadBase) match {
	  case Some(base) => base
	  case None => Base.N // TODO: add better exception handling code - this case shouldn't happen unless there is deletion evidence in the rod
	}

	((readBase != v.getReferenceBase) && v.getRangeLength == 0)
      })

      val productMismatch = if (mismatchBases.length != 0) {
        mismatchBases.map ((base) => {
          val epsilon = Phred.phredToProbability ((base.getMapQuality + base.getSangerQuality) / 2)
          
	  pow ((ploidy - g) * (1 - epsilon) + g * epsilon, base.getCountAtPosition.toDouble)
        }).reduce (_ * _)
      } else {
        0.0
      }

      likelihood (g) = productMatch * productMismatch / pow (ploidy.toDouble, k.toDouble)
    }
    
    likelihood
  }

  /**
   * Calls a SNP for a single pileup, if there is sufficient evidence of a mismatch from the reference.
   * Takes in a single pileup rod from a single sample at a single locus. For simplicity, we assume that
   * all sites are biallelic.
   * 
   * @param[in] pileup List of pileups. Should only contain one rod.
   * @return List of variants seen at site. List can contain 0 or 1 elements - value goes to flatMap.
   */
  protected def callSNP (pileup: List[ADAMPileup]): (List[ADAMVariant], List[ADAMGenotype]) = {
	
    val loci = pileup.head.getPosition
    log.info ("Calling pileup at " + loci)
    
    // get a count of the total number of bases that are a mismatch with the reference
    val nonRefBaseCount = pileup.filter (r => r.getReadBase != r.getReferenceBase)
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
    val maxNonRefBase = if (!nonRefBaseCount.isEmpty) {
      nonRefBaseCount.reduce (pickMaxBase)._1
    } else {
      Base.N // TODO: add better exception handling code
    }

    // score off of rod info
    var likelihood = scoreGenotypeLikelihoods (pileup) 
      
    // compensate likelihoods by minor allele frequency - eqn 19 from source
    likelihood (0) = likelihood (0) * pow ((1.0 - minorAlleleFrequency), 2.0)
    likelihood (1) = likelihood (1) * minorAlleleFrequency * (1.0 - minorAlleleFrequency)
    likelihood (2) = likelihood (2) * pow (minorAlleleFrequency, 2.0)

    // write calls to list
    writeCallInfo (pileup.head, likelihood, maxNonRefBase)
  }

  /**
   * Call variants using simple pileup based SNP calling algorithm.
   *
   * @param[in] pileupGroups An RDD containing lists of pileups.
   * @return An RDD containing called variants.
   */
  override def call (pileups: RDD [ADAMRod]): RDD [(ADAMVariant, List[ADAMGenotype])] = {

    log.info (pileups.count.toString + " rods to call.")

    log.info ("Calling SNPs on pileups and flattening.")
    pileups.map (_.pileups)
      .map (callSNP)
      .filter (_._1.length == 1)
      .map (kv => (kv._1 (0), kv._2))
  }

  override def isCallable (): Boolean = true
}


