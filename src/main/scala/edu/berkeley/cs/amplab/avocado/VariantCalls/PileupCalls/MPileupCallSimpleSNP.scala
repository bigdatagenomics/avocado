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
import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup, Base, ADAMGenotype, VariantType}
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
 * for multiple samples. At the current point in time, we assume that we are running on a diploid organism,
 * and that all sites are biallelic.
 */
class MPileupCallSimpleSNP extends PileupCallSimpleSNP {
  
  override val callName = "MultipleSamplesSimpleSNP"

  /**
   * Calls a SNP for a single pileup, if there is sufficient evidence of a mismatch from the reference.
   * Takes in a single pileup rod from a single sample at a single locus. For simplicity, we assume that
   * all sites are biallelic.
   * 
   * @param[in] pileup List of pileups. Should only contain one rod.
   * @return List of variants seen at site. List can contain 0 or 1 elements - value goes to flatMap.
   */
  protected def callSNP (pileup: ADAMRod): List[ADAMGenotype] = {
	
    val loci = pileup.position
    log.info ("Calling pileup at " + loci)
    
    // get a count of the total number of bases that are a mismatch with the reference
    val nonRefBaseCount = pileup.pileups.filter (r => r.getReadBase != r.getReferenceBase)
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

    val samples = pileup.splitBySamples

    // score off of rod info
    val likelihoods: Array[Array[(Double, Double, Double)]] = samples.map(_.pileups)
      .map(scoreGenotypeLikelihoods)
      .map(v => Array((v(0), v(1), v(2))))
      .toArray
      
    // run maf EM
    val majAlleleFrequency = EMForAlleles.emForMAF(Array(1.0 - minorAlleleFrequency),
                                      likelihoods)

    def compensate (likelihood: Array[(Double, Double, Double)], maf: Double): Array[Double] = {
      // compensate likelihoods by major allele frequency - eqn 19 from source
      Array(likelihood (0)._1 * pow ((1.0 - maf), 2.0),
            likelihood (0)._2 * 2.0 * maf * (1.0 - maf),
            likelihood (0)._3 * pow (maf, 2.0))
    }

    // compensate genotypes by maf
    val compensatedLikelihoods = likelihoods.map(compensate(_, majAlleleFrequency (0)))

    // genotype/variant lists
    var g = List[ADAMGenotype]()
    
    // loop over samples and write calls to list
    for (i <- 0 until likelihoods.length) {
      val l = MutableList(compensatedLikelihoods(i):_*)

      val sg: List[ADAMGenotype] = writeCallInfo (samples(i).pileups.head, l, maxNonRefBase)
    
      g = g ::: sg
    }

    g
  }

  /**
   * Call variants using simple pileup based SNP calling algorithm.
   *
   * @param[in] pileupGroups An RDD containing lists of pileups.
   * @return An RDD containing called variants.
   */
  override def call (pileups: RDD [ADAMRod]): RDD [ADAMGenotype] = {

    log.info (pileups.count.toString + " rods to call.")

    log.info ("Calling SNPs on pileups and flattening.")
    pileups.map (callSNP)
      .filter (_.length != 0)
      .flatMap ((p: List[ADAMGenotype]) => p)
  }

  override def isCallable (): Boolean = true
}


