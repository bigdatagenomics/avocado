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

import breeze.stats.distributions.Binomial
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.adam.rdd.variant.{
  GenotypeRDD,
  VariantContextRDD
}
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.util.LogUtils
import org.bdgenomics.formats.avro.{
  Genotype,
  GenotypeAllele,
  Variant,
  VariantAnnotation,
  VariantCallingAnnotations
}
import scala.collection.JavaConversions._
import scala.math.log

/**
 * Jointly calls variants and computes variant annotations.
 */
object JointAnnotatorCaller extends Serializable {

  private val MINUS_TEN_DIV_LOG10 = -10.0 / log(10.0)

  /**
   * Jointly calls variants and computes variant annotations.
   *
   * @param genotypes The genotypes to jointly process.
   * @return Returns a squared off and annotated set of variant contexts.
   */
  def apply(genotypes: GenotypeRDD): VariantContextRDD = {
    apply(genotypes.toVariantContexts)
  }

  /**
   * Jointly calls variants and computes variant annotations.
   *
   * @param variantContexts The squared off sites to process.
   * @return Returns a squared off and annotated set of variant contexts.
   */
  def apply(variantContexts: VariantContextRDD): VariantContextRDD = {
    variantContexts.transform(_.flatMap(annotateSite))
  }

  /**
   * Jointly calls and annotates a single site.
   *
   * If the site is not variant, we discard the site.
   *
   * @param site The variant site to annotate.
   * @return Returns an annotated and jointly called site.
   */
  private[genotyping] def annotateSite(
    site: VariantContext): Option[VariantContext] = {

    // calculate allele frequency
    // if 0, then we can discard this site
    val minorAlleleFrequency = calculateMinorAlleleFrequency(site)

    if (minorAlleleFrequency <= 0.0) {
      None
    } else {

      // if we have more than one sample, roll up annotations
      val optAnnotations = if (site.genotypes.size > 1) {
        Some(calculateAnnotations(site))
      } else {
        None
      }

      // update the genotypes with our allele frequency
      val updatedGenotypes = recallGenotypes(site,
        minorAlleleFrequency)

      // compute the variant quality
      val variantQuality = computeQuality(updatedGenotypes)

      // make the new variant
      val newVariant = optAnnotations.fold(
        Variant.newBuilder(site.variant.variant))(annotations => {
          Variant.newBuilder(site.variant.variant)
            .setAnnotation(annotations)
        }).setQuality(variantQuality)
        .build

      Some(VariantContext(newVariant, updatedGenotypes))
    }
  }

  /**
   * Calculates the frequency of alternate alleles at a site.
   *
   * @param The site to compute allele frequency for.
   * @return The frequency of alternate alleles.
   */
  private[genotyping] def calculateMinorAlleleFrequency(
    site: VariantContext): Double = {

    val calledAlleles = site.genotypes.map(gt => {
      gt.getAlleles.count(_ != GenotypeAllele.NO_CALL)
    }).sum
    val altAlleles = site.genotypes.map(gt => {
      gt.getAlleles.count(_ == GenotypeAllele.ALT)
    }).sum

    altAlleles.toDouble / calledAlleles.toDouble
  }

  /**
   * Computes rolled up statistical annotations.
   *
   * Computes:
   *
   * - Allelic depth
   * - Strand bias components
   *
   * @param site The site to roll up annotations for.
   * @return A new variant annotation.
   */
  private[genotyping] def calculateAnnotations(
    site: VariantContext): VariantAnnotation = {

    site.genotypes.map(VariantSummary(_))
      .reduce(_.merge(_))
      .toAnnotation(site.variant.variant)
  }

  /**
   * Recalls the genotypes given the estimated allele frequency.
   *
   * @param site The site to rescore genotypes for.
   * @param minorAlleleFrequency The computed alternate allele frequency.
   * @return The new genotype calls to emit.
   */
  private[genotyping] def recallGenotypes(
    site: VariantContext,
    minorAlleleFrequency: Double): Iterable[Genotype] = {

    // if the minor allele frequency is well defined,
    // precompute the prior distributions for all copy number we have
    val priorDistributions: Map[Int, Seq[Double]] = if (minorAlleleFrequency >= 1.0 ||
      minorAlleleFrequency <= 0.0) {
      Map.empty
    } else {
      site.genotypes.map(gt => {
        gt.getAlleles.size
      }).toSet
        .map((cn: Int) => {
          val dist = Binomial(cn, minorAlleleFrequency)

          (cn -> (0 to cn).map(c => {
            dist.logProbabilityOf(c)
          }))
        }).toMap
    }

    site.genotypes.map(recallGenotype(_, priorDistributions))
  }

  private def recallGenotype(gt: Genotype,
                             priors: Map[Int, Seq[Double]]): Genotype = {

    val others = gt.getAlleles.count(a => {
      (a == GenotypeAllele.OTHER_ALT ||
        a == GenotypeAllele.NO_CALL)
    })

    // are the priors defined? if not, just normalize
    // also, skip everything if a no-call/other-alt
    if (priors.nonEmpty && others == 0) {

      // get the prior for our copy number
      val prior = priors(gt.getAlleles.size)

      // merge the prior with the likelihoods
      val nonNormalizedPosterior = new Array[Double](
        gt.getGenotypeLikelihoods.size)
      nonNormalizedPosterior.indices
        .foreach(idx => {
          nonNormalizedPosterior(idx) = (prior(idx) +
            gt.getGenotypeLikelihoods.get(idx))
        })

      // normalize the posterior
      val posterior = LogUtils.logNormalize(nonNormalizedPosterior)
        .map(d => d.toFloat: java.lang.Float)

      // javafy the prior
      val jPrior = prior.map(d => d.toFloat: java.lang.Float)

      // recompute state and quality
      val (state, quality) = BiallelicGenotyper.genotypeStateAndQuality(
        nonNormalizedPosterior)

      // update the current annotations, if available
      val vca = Option(gt.getVariantCallingAnnotations).fold(
        VariantCallingAnnotations.newBuilder)(v => {
          VariantCallingAnnotations.newBuilder(v)
        }).setGenotypePosteriors(posterior.toSeq)
        .setGenotypePriors(jPrior)
        .build

      // build the genotype call array
      val alleles = Seq.fill(state)({
        GenotypeAllele.ALT
      }) ++ Seq.fill(gt.getAlleles.size - state)({
        GenotypeAllele.REF
      })

      Genotype.newBuilder(gt)
        .setVariantCallingAnnotations(vca)
        .setNonReferenceLikelihoods(Seq[java.lang.Double]())
        .setGenotypeQuality(quality.toInt)
        .setAlleles(alleles)
        .build
    } else if (others == 0) {

      // normalize the likelihoods
      val posterior = LogUtils.logNormalize(gt.getGenotypeLikelihoods
        .map(d => d: Double)
        .toArray)
        .map(d => d.toFloat: java.lang.Float)

      // update the current annotations, if available
      val vca = Option(gt.getVariantCallingAnnotations).fold(
        VariantCallingAnnotations.newBuilder)(v => {
          VariantCallingAnnotations.newBuilder(v)
        }).setGenotypePosteriors(posterior.toSeq)
        .build

      Genotype.newBuilder(gt)
        .setVariantCallingAnnotations(vca)
        .setNonReferenceLikelihoods(Seq[java.lang.Double]())
        .build
    } else {
      // WAR for bigdatagenomics/ADAM#1776
      Genotype.newBuilder(gt)
        .setNonReferenceLikelihoods(Seq[java.lang.Double]())
        .build
    }
  }

  /**
   * Computes the variant quality at this site.
   *
   * @param genotypes The genotypes called at this site.
   * @return Returns the Phred scaled probability that this site is variant.
   */
  private[genotyping] def computeQuality(
    genotypes: Iterable[Genotype]): Double = {

    MINUS_TEN_DIV_LOG10 * genotypes.flatMap(gt => {
      Option(gt.getVariantCallingAnnotations)
        .flatMap(vca => {
          Option(vca.getGenotypePosteriors)
        }).flatMap(posteriors => {
          posteriors.headOption
        }).map(_.toDouble)
    }).sum
  }
}
