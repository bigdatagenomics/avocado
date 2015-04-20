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
package org.bdgenomics.avocado.genotyping.annotators

import org.bdgenomics.avocado.algorithms.math.LogBinomial
import org.bdgenomics.avocado.algorithms.math.LogToPhred.log2phred
import org.bdgenomics.avocado.models.AlleleObservation
import org.bdgenomics.formats.avro.{ Variant, VariantCallingAnnotations }

object VariantCallingAnnotator {
  def apply(variant: Variant,
            observations: Iterable[AlleleObservation],
            variantQuality: Double,
            likelihoods: Array[Double]): VariantCallingAnnotator = {
    new VariantCallingAnnotator(variant, observations, variantQuality, likelihoods)
  }
}

class VariantCallingAnnotator(variant: Variant, observations: Iterable[AlleleObservation], variantQuality: Double, likelihoods: Array[Double]) {

  val alt = variant.getAlternateAllele
  val ref = variant.getReferenceAllele
  val altObservations = observations.filter(_.allele == alt)
  val refObservations = observations.filter(_.allele == ref)
  val numAlt = altObservations.size
  val numRef = refObservations.size

  val annotationBuilder = VariantCallingAnnotations.newBuilder()

  def build(): VariantCallingAnnotations = {
    // get optional rank sum
    mapQRankSum().foreach(annotationBuilder.setMqRankSum(_))

    // other values are not optional
    annotationBuilder.setVariantQualityByDepth(variantQualityByDepth())
      .setReadDepth(readDepth())
      .setFisherStrandBiasPValue(fisherStrandBiasValue())
      .setRmsMapQ(rmsMapQ())
      .build()
  }

  def variantQualityByDepth(): Float = {
    (variantQuality / numAlt).toFloat
  }

  def readDepth(): Integer = {
    numAlt + numRef
  }

  // Calculate RankSum as according to https://www.broadinstitute.org/gatk/guide/article?id=473
  def mapQRankSum(): Option[Float] = {

    if (altObservations.isEmpty || refObservations.isEmpty) {
      None
    } else {
      val sortedObservations: Array[AlleleObservation] = (altObservations ++ refObservations).toArray.sortBy(_.mapq)
      //Make (observation, rank) tuples
      val rankedObservations: Array[(AlleleObservation, Float)] = sortedObservations.map((x: AlleleObservation) => (x, 0f))

      for (idx <- 0 until sortedObservations.length) {
        val curMapQ = sortedObservations(idx).mapq

        // Scan above and under for same mapq
        var minIdx = idx
        while (minIdx > 0 && sortedObservations(minIdx - 1).mapq == curMapQ) {
          minIdx = minIdx - 1
        }
        var maxIdx = idx
        while (maxIdx < sortedObservations.size - 1 && sortedObservations(maxIdx + 1).mapq == curMapQ) {
          maxIdx = maxIdx + 1
        }

        val rank = (minIdx + maxIdx) / 2f + 1
        rankedObservations(idx) = rankedObservations(idx).copy(_2 = rank)
      }

      val rankRef = rankedObservations.filter(_._1.allele == ref).map(_._2).reduce(_ + _)
      val rankAlt = rankedObservations.filter(_._1.allele == alt).map(_._2).reduce(_ + _)

      val uRef = rankRef - numRef * (numRef + 1) / 2f
      val uAlt = rankAlt - numAlt * (numAlt + 1) / 2f

      val mu = (numRef * numAlt) / 2f
      val std = math.sqrt((numRef * numAlt * (numRef + numAlt + 1)) / 12f)
      Some(((uAlt - mu) / std).toFloat)
    }
  }

  def fisherStrandBiasValue(): Float = {
    val negativeAlts = altObservations.filter(_.onNegativeStrand).size
    val positiveAlts = altObservations.filter(!_.onNegativeStrand).size
    val negativeRefs = refObservations.filter(_.onNegativeStrand).size
    val positiveRefs = refObservations.filter(!_.onNegativeStrand).size
    val n = positiveAlts + positiveRefs + negativeAlts + negativeRefs
    val logP = LogBinomial.logBinomial(positiveRefs + positiveAlts, positiveRefs) +
      LogBinomial.logBinomial(negativeRefs + negativeAlts, negativeRefs) -
      LogBinomial.logBinomial(n, positiveRefs + negativeRefs)
    log2phred(logP).toFloat
  }

  def rmsMapQ(): Float = {
    val mapqObservations = observations.filter(!_.mapq.isEmpty)
    val meanSquareMapQ = (mapqObservations.map(x => math.pow(x.mapq.get.toDouble, 2)).reduce(_ + _)) / mapqObservations.size.toFloat
    math.sqrt(meanSquareMapQ).toFloat
  }
}
