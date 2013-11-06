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
import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup,ADAMVariant,Base}
import edu.berkeley.cs.amplab.avocado.utils.Phred
import scala.math.pow
import scala.collection.mutable.MutableList

/**
 * Trait for filtering pileups. 
 */
class PileupCallSimpleSNP extends PileupCall ("") {
  
  // get ploidy - if no ploidy, assume human diploid
  val ploidy = 2

  /**
   * Calls a SNP for a single pileup, if mismatch from the reference.
   *
   * @param[in] pileup List of pileups. Should only contain one.
   * @return List of variants seen at site. Can range from 0 to 2.
   */
  def callSNP (pileup: List[ADAMPileup]): List[ADAMVariant] = {
	      
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
    var likelihood = MutableList[Double](3)
    val k = pileup.map (_.getCountAtPosition).reduce (_ + _)

    for (g <- 0 to 2) {
      val productMatch = pileup.filter(v => v.getReadBase == v.getReferenceBase).map ((base) => {
        val epsilon = Phred.phredToProbability ((base.getMapQuality + base.getSangerQuality) / 2)
        
	pow ((ploidy - g) * epsilon + g * (1 - epsilon), base.getCountAtPosition.toDouble)
      }).reduce (_ * _)
     
      val productMismatch = pileup.filter(v => {
	val readBase = Option (v.getReadBase) match {
	  case Some(base) => base
	  case None => Base.N // FIXME
	}

	((readBase != v.getReferenceBase) && v.getRangeLength == 0)
      }).map ((base) => {
        val epsilon = Phred.phredToProbability ((base.getMapQuality + base.getSangerQuality) / 2)
        
	pow ((ploidy - g) * (1 - epsilon) + g * epsilon, base.getCountAtPosition.toDouble)
      }).reduce (_ * _)
      
      likelihood (g) = productMatch * productMismatch / pow (ploidy.toDouble, k.toDouble)
    }

    // simplifying assumption - snps are biallelic
    // will address later
    if (likelihood.indexOf (likelihood.max) == 0) {
      // if max likelihood is homozygous reference, return no variants
      List[ADAMVariant]()
    } else if (likelihood.indexOf (likelihood.max) == 1) {
      // if max likelihood is in position 1, we are heterozygous, return no variants
      List[ADAMVariant]()
    } else {
      // if neither of the above are true, we are homozygous but not reference
      // then, need to find max base pair probability
      List[ADAMVariant]()
    }
      
  }

  /**
   * 
   *
   * @param[in] pileupGroups An RDD containing lists of pileups.
   * @return An RDD containing called variants.
   */
  override def call (pileups: RDD [ADAMPileup]): RDD [ADAMVariant] = {
    return pileups.groupBy (_.getPosition).map (_._2.toList).flatMap (callSNP)
  }
}


