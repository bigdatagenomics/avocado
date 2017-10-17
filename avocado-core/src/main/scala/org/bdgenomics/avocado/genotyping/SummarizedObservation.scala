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
package org.bdgenomics.avocado.genotyping

import org.bdgenomics.avocado.models.Observation

/**
 * Represents a scored allele, before statistics are calculated.
 *
 * This representation is used with the ScoredObservation class/object,
 * which is used to precompute likelihoods and other stats, which saves
 * on runtime.
 *
 * @param isRef Was this observation of a reference allele?
 * @param forwardStrand Was this on the forward strand?
 * @param optQuality What was the quality score of the bases in the observation,
 *   if known.
 * @param mapQ What was the mapping quality of the bases in the observation?
 * @param isOther Was this observation of an other-alt allele?
 * @param isOther Was this observation of an unknown allele?
 */
private[genotyping] case class SummarizedObservation(isRef: Boolean,
                                                     forwardStrand: Boolean,
                                                     optQuality: Option[Int],
                                                     mapQ: Int,
                                                     copyNumber: Option[Int] = None,
                                                     isOther: Boolean = false,
                                                     isNonRef: Boolean = false) {

  /**
   * Associates a copy number with this observation.
   *
   * @param copyNumber The copy number of this site.
   * @return Returns a new observation with known copy number.
   */
  def addCopyNumber(copyNumber: Int): SummarizedObservation = {
    copy(copyNumber = Some(copyNumber))
  }

  /**
   * @return Returns an observation that agrees with the reference allele.
   */
  def asRef: SummarizedObservation = {
    SummarizedObservation(true,
      forwardStrand,
      optQuality,
      mapQ)
  }

  /**
   * @return Returns an observation that agrees with the alternate allele.
   */
  def asAlt: SummarizedObservation = {
    SummarizedObservation(false,
      forwardStrand,
      optQuality,
      mapQ)
  }

  /**
   * @return Returns an observation that is another alternate allele.
   */
  def otherAlt: SummarizedObservation = {
    SummarizedObservation(false,
      forwardStrand,
      optQuality,
      mapQ,
      isOther = true)
  }

  /**
   * @return Returns an observation that is an unknown allele.
   */
  def nullOut: SummarizedObservation = {
    SummarizedObservation(false,
      forwardStrand,
      optQuality,
      mapQ,
      isNonRef = true)
  }

  /**
   * @param scores The previously scored observations.
   * @return Returns a fully scored observation.
   */
  def toObservation(scores: Seq[ScoredObservation]): Observation = {
    val optScore = scores.find(score => {
      score.isRef == isRef &&
        score.forwardStrand == forwardStrand &&
        score.optQuality == optQuality &&
        score.mapQ == mapQ &&
        score.isOther == isOther &&
        score.isNonRef == isNonRef &&
        copyNumber.fold(true)(cn => cn == score.copyNumber)
    })
    assert(optScore.isDefined)
    val score = optScore.get

    Observation(score.alleleForwardStrand,
      score.otherForwardStrand,
      score.squareMapQ,
      score.referenceLogLikelihoods,
      score.alleleLogLikelihoods,
      score.otherLogLikelihoods,
      score.nonRefLogLikelihoods,
      score.alleleCoverage,
      score.otherCoverage,
      score.totalCoverage,
      score.isRef,
      score.copyNumber)
  }
}
