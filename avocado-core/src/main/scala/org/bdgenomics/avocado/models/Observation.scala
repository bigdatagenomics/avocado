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
package org.bdgenomics.avocado.models

/**
 * A generic class that stores likelihoods and simple annotations.
 *
 * @param forwardStrand The number of reads observed on the forward strand.
 * @param squareMapQ The sum of the squares of the mapping qualities observed.
 * @param alleleLogLikelihoods The log likelihoods that 0...n copies of this allele were
 *   observed.
 * @param otherLogLikelihoods The log likelihoods that 0...n copies of another allele were
 *   observed.
 * @param coverage The total number of reads observed.
 */
case class Observation(forwardStrand: Int,
                       squareMapQ: Double,
                       alleleLogLikelihoods: Array[Double],
                       otherLogLikelihoods: Array[Double],
                       coverage: Int) {

  protected val copyNumber = alleleLogLikelihoods.length
  assert(copyNumber == otherLogLikelihoods.length &&
    copyNumber > 0)
  assert(squareMapQ >= 0.0)
  assert(coverage > 0)
  assert(forwardStrand >= 0)

  /**
   * Merges two observations.
   *
   * @note This method destructively updates the first observation by modifying
   *   the underlying likelihood arrays in place.
   *
   * @param obs Observation to merge with.
   * @return Returns a new observation that is the sum of the two input
   *   observations.
   */
  def merge(obs: Observation): Observation = {
    assert(copyNumber == obs.copyNumber)

    (0 until copyNumber).foreach(i => {
      alleleLogLikelihoods(i) += obs.alleleLogLikelihoods(i)
      otherLogLikelihoods(i) += obs.otherLogLikelihoods(i)
    })

    Observation(forwardStrand + obs.forwardStrand,
      squareMapQ + obs.squareMapQ,
      alleleLogLikelihoods,
      otherLogLikelihoods,
      coverage + obs.coverage)
  }
}
