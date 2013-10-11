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

import spark.{RDD,SparkContext}
import edu.berkeley.cs.amplab.adam.util.{Pileup,PileupTraversable}
import edu.berkeley.cs.amplab.adam.avro.ADAMVariant
import edu.berkeley.cs.amplab.avocado.utils.Phred

/**
 * Trait for filtering pileups. 
 */
class PileupCallSimpleSNP extends PileupCall {
  
  // get ploidy - if no ploidy, assume human diploid
  try {
    val ploidy = config [Int]("ploidy")
  } catch {
    case NoSuchElementException nse: val ploidy = 2
  }

  /**
   * Calls a SNP for a single pileup, if mismatch from the reference.
   *
   * @param[in] pileup List of pileups. Should only contain one.
   * @return List of variants seen at site. Can range from 0 to 2.
   */
  def callSNP (pileup: List[Pileup]): List[ADAMVariant] = {
    require (pileup.length = 1, 
	     {println ("Can only provide one pileup to SNP caller. Pileup provided covers %d to %d.", 
		       pileup.reduce (min (_.referencePosition, _.referencePosition)),
		       pileup.reduce (max (_.referencePosition, _.referencePosition)))})
    require (pileup.head.deletes.length == 0 &&
	     pileup.head.insertions.length = 0, 
	     {println ("Can only call SNP on pileups without INDEL evidence. Pileup at %d.",
		       pileup.head.referencePosition)}
	      
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
    var likelihood = Array[Int](3)
    val k = pileup.head.matches.length + pileup.head.mismatches.length

    for (g <- 0 to 2) {
      val productMatch = pileup.head.matches.map ((base) => {
	val epsilon = phredToProbability ((base.mapQ + base.qual) / 2.0)
	return (ploidy - g) * epsilon + g * (1 - epsilon)
      }).reduce (_ * _)
      val productMismatch = pileup.head.mismatches.map ((base) => {
	val epsilon = phredToProbability ((base.mapQ + base.qual) / 2.0)
	return (ploidy - g) * (1 - epsilon) + g * epsilon
      }).reduce (_ * _)
      
      likelihood (g) = productMatch * productMismatch / pow (ploidy, k)
    }

    if (likelihood.indexOf (likelihood.max) == 0) {
      return List[ADAMVariant](0) // if max likelihood is homozygous reference, return no variants
    }
      
  }

  /**
   * 
   *
   * @param[in] pileupGroups An RDD containing lists of pileups.
   * @return An RDD containing called variants.
   */
  override def call (pileupGroups: RDD [(void, List[Pileup])]): RDD [ADAMVariant] = {
    return pileupGroups.flatMap (callSNP)
  }
}

