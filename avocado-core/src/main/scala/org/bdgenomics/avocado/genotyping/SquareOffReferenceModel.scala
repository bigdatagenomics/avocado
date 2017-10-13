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

import org.bdgenomics.adam.models.{ ReferenceRegion, VariantContext }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.variant.{
  GenotypeRDD,
  VariantRDD,
  VariantContextRDD
}
import org.bdgenomics.formats.avro.{ Genotype, GenotypeAllele, Variant }
import scala.annotation.tailrec

/**
 * Squares off a set of genotypes with reference models.
 *
 * Many joint genotyping workflows use a "Genome VCF" (gVCF) based approach to
 * incrementally compute genotype likelihoods across their dataset. In this
 * methodology, we generate genotype likelihoods at all positions in all
 * samples. For sites where we do not see evidence of a variant, we compute a
 * "reference model", which is a set of genotype likelihoods assuming that we
 * saw an unknown alternate allele. These likelihoods are then used in a joint
 * genotyping step.
 *
 * The alternative to this approach is to discover variants across all samples
 * simultaneously, and to then score these variants. This approach is generally
 * considered too computationally expensive for large cohorts.
 *
 * This singleton object "squares off" the reference model by discovering all
 * sites where we called a variant in at least one sample, joining these
 * discovered variants back against the input genotypes, and then excising the
 * genotype likelihoods from the reference models.
 */
object SquareOffReferenceModel {

  /**
   * Squares off genotypes containing both called sites and reference models.
   *
   * @param genotypes Genotypes containing both called sites and reference models.
   * @return A set of variant contexts where at least one copy of the alternate
   *   allele was called across all samples, with genotype likelihood models for
   *   all samples that had data at the site.
   */
  def apply(genotypes: GenotypeRDD): VariantContextRDD = {

    val variants = extractVariants(genotypes)

    // join variants back against genotypes
    val sites = variants.shuffleRegionJoinAndGroupByLeft(genotypes)

    sites.transmute(_.map(s => squareOffSite(s._1, s._2)))
  }

  /**
   * Trims the rightmost end of a ref/alt allele pair.
   *
   * Needed to canonicalize alleles across samples. If the two alleles in a
   * ref/alt allele pair end with the same suffix, trims those bases.
   *
   * @param ref The reference allele.
   * @param alt The alternate allele.
   * @return Returns the number of bases to trim from the two alleles, if any.
   */
  private[genotyping] def trimRight(ref: String,
                                    alt: String): Int = {

    @tailrec def trim(refIter: Iterator[Char],
                      altIter: Iterator[Char],
                      trimmed: Int = 0): Int = {
      if (refIter.next != altIter.next) {
        trimmed
      } else if (!(refIter.hasNext && altIter.hasNext)) {
        trimmed
      } else {
        trim(refIter, altIter, trimmed + 1)
      }
    }

    if (ref.isEmpty || alt.isEmpty) {
      0
    } else {
      trim(ref.reverse.toIterator,
        alt.reverse.toIterator)
    }
  }

  /**
   * Discovers variant sites from the reference model genotypes.
   *
   * @param genotypes Genotypes containing both called sites and reference models.
   * @return Returns sites where a variant was seen in at least one sample.
   */
  private[genotyping] def extractVariants(genotypes: GenotypeRDD): VariantRDD = {

    genotypes.transmute(_.flatMap(gt => {

      if (gt.getAlleles.contains(GenotypeAllele.ALT) &&
        gt.getVariant.getAlternateAllele != null) {

        // trim the alleles so that we have a canonical result
        val toTrimFromRight = trimRight(gt.getVariant.getReferenceAllele,
          gt.getVariant.getAlternateAllele)

        Some(Variant.newBuilder
          .setContigName(gt.getContigName)
          .setStart(gt.getStart)
          .setEnd(gt.getEnd - toTrimFromRight)
          .setReferenceAllele(gt.getVariant
            .getReferenceAllele
            .dropRight(toTrimFromRight))
          .setAlternateAllele(gt.getVariant
            .getAlternateAllele
            .dropRight(toTrimFromRight))
          .build)
      } else {
        None
      }
    })).transformDataset(_.dropDuplicates("start",
      "end",
      "contigName",
      "referenceAllele",
      "alternateAllele"))
  }

  /**
   * Processes all of the genotypes overlapping a variant.
   *
   * @param variant The variant we are interested in.
   * @param genotypes The genotypes overlapping this variant, including both
   *   called genotypes and reference model genotypes.
   * @return Returns genotypes whose likelihoods model the variant we are
   *   interested in.
   */
  private[genotyping] def squareOffSite(
    variant: Variant,
    genotypes: Iterable[Genotype]): VariantContext = {

    VariantContext(
      variant,
      genotypes.groupBy(_.getSampleId)
        .flatMap(kv => {
          val (_, sampleGenotypes) = kv

          hasScoredVariant(variant, sampleGenotypes)
            .orElse(exciseGenotypeFromReferenceModel(variant, sampleGenotypes))
        }))
  }

  /**
   * Picks the genotype that contains the scored variant, if one exists.
   *
   * @param variant The variant we are interested in.
   * @param genotypes The genotypes we saw from a single sample that covered the
   *   variant of interest.
   * @return If one exists, the genotype corresponding to the variant we are
   *   interested in.
   */
  private[genotyping] def hasScoredVariant(
    variant: Variant,
    genotypes: Iterable[Genotype]): Option[Genotype] = {

    genotypes.find(gt => {
      gt.getStart == variant.getStart &&
        gt.getEnd == variant.getEnd &&
        gt.getContigName == variant.getContigName &&
        gt.getVariant.getReferenceAllele == variant.getReferenceAllele &&
        gt.getVariant.getAlternateAllele == variant.getAlternateAllele
    })
  }

  /**
   * Given reference models, extracts the likelihood model for a variant.
   *
   * @param variant The variant we are interested in.
   * @param genotypes The genotypes we saw from a single sample that covered the
   *   variant of interest.
   * @return The likelihood model for the variant, if it can be computed.
   */
  private[genotyping] def exciseGenotypeFromReferenceModel(
    variant: Variant,
    genotypes: Iterable[Genotype]): Option[Genotype] = {

    val rr = ReferenceRegion(variant)

    genotypes.find(gt => {
      val gtRr = ReferenceRegion(gt)

      gtRr.overlaps(rr)
    }).map(gt => {
      Genotype.newBuilder(gt)
        .setStart(variant.getStart)
        .setEnd(variant.getEnd)
        .setVariant(variant)
        .setGenotypeLikelihoods(gt.getNonReferenceLikelihoods)
        .build
    })
  }
}
