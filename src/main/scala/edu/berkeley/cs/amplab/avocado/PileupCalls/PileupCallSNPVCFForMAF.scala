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
import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup,ADAMVariant,Base,ADAMGenotype,VariantType}
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
class PileupCallSNPVCFForMAF extends PileupCallSimpleSNP {
  
  override val callName = "SNPVCFForMAF"

  /**
   * Calls a SNP for a single pileup, if there is sufficient evidence of a mismatch from the reference.
   * Takes in a single pileup rod from a single sample at a single locus. For simplicity, we assume that
   * all sites are biallelic.
   * 
   * @param[in] pileup List of pileups. Should only contain one rod.
   * @return List of variants seen at site. List can contain 0 or 1 elements - value goes to flatMap.
   */
  override protected def callSNP (pileup: List[ADAMPileup]): (List[ADAMVariant], List[ADAMGenotype]) = {
	
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
    val maxNonRefBase = nonRefBaseCount.reduce (pickMaxBase)._1

    // score off of rod info
    var likelihood = scoreGenotypeLikelihoods (pileup) 
      
    // compensate likelihoods by minor allele frequency - get from vcf
    likelihood (0) = likelihood (0) * pow ((1.0), 2.0)
    likelihood (1) = likelihood (1) * (1.0)
    likelihood (2) = likelihood (2) * pow (1.0, 2.0)

    // write calls to list
    writeCallInfo (pileup.head, likelihood, maxNonRefBase)
  }
}


