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

import org.bdgenomics.adam.avro.{ Base, ADAMGenotype, ADAMPileup }
import org.bdgenomics.adam.models.{ ADAMRod, ADAMVariantContext }
import org.bdgenomics.avocado.calls.VariantCallCompanion
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import scala.collection.mutable.MutableList
import scala.collection.JavaConversions._
import scala.math.pow

object MPileupCallSimpleSNP extends VariantCallCompanion {

  val callName = "MultipleSamplesSimpleSNP"

  def apply(stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): PileupCallSimpleSNP = {

    val ploidy = config.getInt("ploidy", 2)
    val minorAlleleFrequency = config.getDouble("minorAlleleFrequency", 0.1)

    new MPileupCallSimpleSNP(ploidy, minorAlleleFrequency)
  }
}

/**
 * Class to call SNPs. Implements the SNP genotyping methods described in:
 *
 * Li, Heng. "A statistical framework for SNP calling, mutation discovery, association mapping and population
 * genetical parameter estimation from sequencing data." Bioinformatics 27.21 (2011): 2987-2993.
 *
 * for multiple samples. At the current point in time, we assume that we are running on a diploid organism,
 * and that all sites are biallelic.
 */
class MPileupCallSimpleSNP(ploidy: Int,
                           minorAlleleFrequency: Double) extends PileupCallSimpleSNP(ploidy) {

  override val companion = MPileupCallSimpleSNP

  /**
   * Calls a SNP for a single pileup, if there is sufficient evidence of a mismatch from the reference.
   * Takes in a single pileup rod from a single sample at a single locus. For simplicity, we assume that
   * all sites are biallelic.
   *
   * @param[in] pileup List of pileups. Should only contain one rod.
   * @return List of variants seen at site. List can contain 0 or 1 elements - value goes to flatMap.
   */
  protected def callSNP(pileup: ADAMRod): List[ADAMVariantContext] = {
    val samples = pileup.splitBySamples

    // score off of rod info
    val likelihoods: Array[Array[(Double, Double, Double)]] = samples.map(_.pileups)
      .map(scoreGenotypeLikelihoods)
      .map(v => Array((v(0), v(1), v(2))))
      .toArray

    // run maf EM
    val majAlleleFrequency = EMForAlleles.emForMAF(Array(1.0 - minorAlleleFrequency),
      likelihoods)

    def compensate(likelihood: Array[(Double, Double, Double)], maf: Double): List[Double] = {
      // compensate likelihoods by major allele frequency - eqn 19 from source
      List(likelihood(0)._1 * pow(maf, 2.0),
        likelihood(0)._2 * 2.0 * maf * (1.0 - maf),
        likelihood(0)._3 * pow((1.0 - maf), 2.0))
    }

    // compensate genotypes by maf
    val compensatedLikelihoods = likelihoods.map(compensate(_, majAlleleFrequency(0)))

    // genotype/variant lists
    var g = List[ADAMGenotype]()

    // get most often seen non-reference base
    val maxNonRefBase = getMaxNonRefBase(pileup.pileups)

    // loop over samples and write calls to list
    for (i <- 0 until likelihoods.length) {
      val l = compensatedLikelihoods(i)

      val sg: List[ADAMGenotype] = writeCallInfo(samples(i).pileups.head, l, maxNonRefBase)

      g = g ::: sg
    }

    genotypesToVariantContext(g, samples.length)
  }

  override def isCallable(): Boolean = true
}

