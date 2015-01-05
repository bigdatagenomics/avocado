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
package org.bdgenomics.avocado.algorithms.em

import org.bdgenomics.avocado.algorithms.math._
import scala.math.{ abs, pow, log }

object EMForAlleles {

  /**
   * Main MAF EM function
   *  IN: Phi - an initial MAF vector of length number of SNPs
   *      GL - Array of arrays of likelihood triples P( D | g )
   *          (note these are NOT multiplied by P(g | phi)! )
   *  OUT: Phi - ML estimate of MAF's across SNPs
   *
   *  Note: GL is currently an array of (numSnps) arrays of length (numInds),
   */
  def emForMAF(genotypeLikelihoods: Array[Array[Double]],
               logReferenceFrequencyEstimate: Double,
               maxIterations: Option[Int] = None,
               targetTolerance: Option[Double] = None): Double = {
    assert(maxIterations.isDefined || targetTolerance.isDefined,
      "Must define at least one of the iteration or tolerance limits.")

    // map target tolerance to log
    val logTargetTolerance = targetTolerance.map(log(_))

    // M is the total number of chromosomes in all samples
    val samplePloidy = genotypeLikelihoods.map(_.length - 1)
    val logM = log(samplePloidy.sum.toDouble)
    val ploidies = samplePloidy.distinct.toArray

    // loop until convergence
    var logPsi = logReferenceFrequencyEstimate
    var lastLogPsi = logPsi
    var iter = 0

    // precompute log of genotype state
    val gtStateMap = ploidies.map(g => {
      (g + 1, (0 to g).map(i => log(i.toDouble)).toArray)
    }).toMap

    do {
      // carry over psi from previous iteration
      lastLogPsi = logPsi

      // calculate the new prior distributions per ploidy
      val ploidyDistributionMap = ploidies.map(m => {
        (m + 1, LogBinomial.calculateLogProbabilities(logPsi, m))
      }).toMap

      // per sample, calculate the contribution of each genotype state
      // then, sum these contributions together and normalize
      logPsi = LogUtils.sumLogProbabilities(genotypeLikelihoods.map(gls => {
        // the length of the genotype likelihood array is equal to the sample ploidy plus 1
        val ploidyP1 = gls.length

        // from this, recover the state prior probabilities and log of genotype state
        val prior = ploidyDistributionMap(ploidyP1)
        val logGt = gtStateMap(ploidyP1)

        // loop to sum
        var num = 0.0
        var denom = 0.0
        (0 until ploidyP1).foreach(i => {
          val contribution = gls(i) + prior(i)

          if (i == 0) {
            num = logGt(i) + contribution
            denom = contribution
          } else {
            num = LogUtils.logSum(num, logGt(i) + contribution)
            denom = LogUtils.logSum(denom, contribution)
          }
        })

        num - denom
      })) - logM

      // increment iteration count
      iter += 1
    } while (logTargetTolerance.fold(true)(_ < abs(logPsi - lastLogPsi)) &&
      maxIterations.fold(true)(_ > iter) &&
      logPsi < 0.0)

    logPsi
  }
}
