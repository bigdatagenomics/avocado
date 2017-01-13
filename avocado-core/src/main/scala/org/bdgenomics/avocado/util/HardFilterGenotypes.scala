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
package org.bdgenomics.avocado.util

import org.bdgenomics.adam.rdd.variant.GenotypeRDD
import org.bdgenomics.formats.avro.{
  Genotype,
  GenotypeAllele,
  VariantCallingAnnotations
}
import scala.collection.JavaConversions._

private[avocado] trait HardFilterGenotypesArgs extends Serializable {

  /**
   * The minimum genotype quality to emit.
   */
  var minQuality: Int

  /**
   * The minimum quality by depth to allow for heterozygous SNPs.
   *
   * Set to a negative value to omit.
   */
  var minHetSnpQualityByDepth: Float

  /**
   * The minimum quality by depth to allow for homozygous SNPs.
   *
   * Set to a negative value to omit.
   */
  var minHomSnpQualityByDepth: Float

  /**
   * The minimum quality by depth to allow for heterozygous INDELs.
   *
   * Set to a negative value to omit.
   */
  var minHetIndelQualityByDepth: Float

  /**
   * The minimum quality by depth to allow for homozygous INDELs.
   *
   * Set to a negative value to omit.
   */
  var minHomIndelQualityByDepth: Float

  /**
   * The max Phred-scaled strand bias p-value to allow for SNPs.
   *
   * Set to a negative value to omit.
   */
  var maxSnpPhredStrandBias: Float

  /**
   * The minimum coverage to allow for SNPs.
   *
   * Set to a negative value to omit.
   */
  var minSnpDepth: Int

  /**
   * The max coverage to allow for SNPs.
   *
   * Set to a negative value to omit.
   */
  var maxSnpDepth: Int

  /**
   * The max Phred-scaled strand bias p-value to allow for INDELs.
   *
   * Set to a negative value to omit.
   */
  var maxIndelPhredStrandBias: Float

  /**
   * The minimum RMS mapping quality to allow for SNPs.
   *
   * Set to a negative value to omit.
   */
  var minSnpRMSMappingQuality: Float

  /**
   * The minimum RMS mapping quality to allow for INDELs.
   *
   * Set to a negative value to omit.
   */
  var minIndelRMSMappingQuality: Float

  /**
   * The minimum coverage to allow for INDELs.
   *
   * Set to a negative value to omit.
   */
  var minIndelDepth: Int

  /**
   * The max coverage to allow for INDELs.
   *
   * Set to a negative value to omit.
   */
  var maxIndelDepth: Int
}

/**
 * Reifies genotypes down to confident alt calls.
 */
private[avocado] object HardFilterGenotypes extends Serializable {

  /**
   * Applies hard filters to a GenotypeRDD.
   *
   * @param grdd GenotypeRDD to filter.
   * @param args The hard filter configuration to apply.
   * @return A new GenotypeRDD of hard filtered genotypes.
   */
  def apply(grdd: GenotypeRDD,
            args: HardFilterGenotypesArgs): GenotypeRDD = {

    // make snp and indel filters
    val snpFilters = buildSnpHardFilters(args)
    val indelFilters = buildIndelHardFilters(args)

    // flat map the filters over the genotype rdd
    val minQuality = args.minQuality
    grdd.transform(rdd => {
      rdd.flatMap(filterGenotype(_,
        minQuality,
        snpFilters,
        indelFilters))
    })
  }

  /**
   * @param args The hard filter configuration to apply.
   * @return A collection of hard filtering functions to apply to SNPs.
   */
  private[util] def buildSnpHardFilters(
    args: HardFilterGenotypesArgs): Iterable[Genotype => Option[String]] = {
    val optHetQdFilter: Option[Genotype => Option[String]] = Option(args.minHetSnpQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, hom = false))
    val optHomQdFilter: Option[Genotype => Option[String]] = Option(args.minHomSnpQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, hom = true))
    val optMqFilter: Option[Genotype => Option[String]] = Option(args.minSnpRMSMappingQuality)
      .filter(_ > 0.0)
      .map(minMq => hardFilterRMSMapQ(_, minMq))
    val optFsFilter: Option[Genotype => Option[String]] = Option(args.maxSnpPhredStrandBias)
      .filter(_ > 0.0)
      .map(minFs => hardFilterStrandBias(_, minFs))
    val optMinDpFilter: Option[Genotype => Option[String]] = Option(args.minSnpDepth)
      .filter(_ > 0)
      .map(minDp => hardFilterMinDepth(_, minDp))
    val optMaxDpFilter: Option[Genotype => Option[String]] = Option(args.maxSnpDepth)
      .filter(_ > 0)
      .map(maxDp => hardFilterMaxDepth(_, maxDp))

    Iterable(optHetQdFilter,
      optHomQdFilter,
      optMqFilter,
      optFsFilter,
      optMinDpFilter,
      optMaxDpFilter).flatten
  }

  /**
   * @param args The hard filter configuration to apply.
   * @return A collection of hard filtering functions to apply to INDELs.
   */
  private[util] def buildIndelHardFilters(
    args: HardFilterGenotypesArgs): Iterable[Genotype => Option[String]] = {
    val optHetQdFilter: Option[Genotype => Option[String]] = Option(args.minHetIndelQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, hom = false))
    val optHomQdFilter: Option[Genotype => Option[String]] = Option(args.minHomIndelQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, hom = true))
    val optMqFilter: Option[Genotype => Option[String]] = Option(args.minIndelRMSMappingQuality)
      .filter(_ > 0.0)
      .map(minMq => hardFilterRMSMapQ(_, minMq))
    val optFsFilter: Option[Genotype => Option[String]] = Option(args.maxIndelPhredStrandBias)
      .filter(_ > 0.0)
      .map(minFs => hardFilterStrandBias(_, minFs))
    val optMinDpFilter: Option[Genotype => Option[String]] = Option(args.minIndelDepth)
      .filter(_ > 0)
      .map(minDp => hardFilterMinDepth(_, minDp))
    val optMaxDpFilter: Option[Genotype => Option[String]] = Option(args.maxIndelDepth)
      .filter(_ > 0)
      .map(maxDp => hardFilterMaxDepth(_, maxDp))

    Iterable(optHetQdFilter,
      optHomQdFilter,
      optMqFilter,
      optFsFilter,
      optMinDpFilter,
      optMaxDpFilter).flatten
  }

  /**
   * @param genotype call to evaluate.
   * @return Returns false for calls that are homozygous ref.
   */
  private[util] def filterRefCalls(genotype: Genotype): Boolean = {
    !genotype.getAlleles.forall(_ == GenotypeAllele.REF)
  }

  /**
   * @param genotype call to evaluate.
   * @param minQuality The minimum genotype quality to emit.
   * @return Returns false for calls that are below the minimum quality.
   */
  private[util] def filterQuality(genotype: Genotype,
                                  minQuality: Int): Boolean = {
    genotype.getGenotypeQuality > minQuality
  }

  /**
   * Filters whether a genotype should be emitted or not.
   *
   * This is a prefilter before the hard filters are applied. If a genotype call
   * is not an alt call or is very low confidence, we shouldn't treat it as an
   * emitted variant. Rather, it is just a site we happened to observe.
   *
   * @param genotype call to evaluate.
   * @param minQuality The minimum genotype quality to emit.
   * @return Returns false for calls that are hom ref and/or low quality.
   */
  def emitGenotypeFilter(genotype: Genotype,
                         minQuality: Int): Boolean = {
    filterRefCalls(genotype) && filterQuality(genotype, minQuality)
  }

  /**
   * @param genotype call to evaluate.
   * @param minQualityByDepth The minimum quality/depth value to accept.
   * @param hom If true, this filter should be applied to homozygous variants.
   * @return If this genotype has a low quality for the depth of the site, we
   *   return an optional filter string. Else, we return a None.
   */
  private[util] def hardFilterQualityByDepth(genotype: Genotype,
                                             minQualityByDepth: Float,
                                             hom: Boolean = false): Option[String] = {
    val gtIsHom = genotype.getAlleles.forall(_ == GenotypeAllele.ALT)
    if ((!gtIsHom && !hom) || (gtIsHom && hom)) {
      Option(genotype.getReadDepth)
        .flatMap(depth => {
          Option(genotype.getGenotypeQuality)
            .map(gq => gq.toFloat / depth.toFloat)
        }).flatMap(qualByDepth => {
          if (qualByDepth < minQualityByDepth) {
            Some("QD<%1.1f".format(minQualityByDepth))
          } else {
            None
          }
        })
    } else {
      None
    }
  }

  /**
   * @param genotype call to evaluate.
   * @param maxPhredStrandBias The maximum Phred strand bias score to allow.
   * @return If this genotype shows significant strand bias, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterStrandBias(genotype: Genotype,
                                         maxPhredStrandBias: Float): Option[String] = {
    Option(genotype.getVariantCallingAnnotations)
      .flatMap(vca => Option(vca.getFisherStrandBiasPValue))
      .flatMap(fs => {
        if (fs > maxPhredStrandBias) {
          Some("FS>%1.1f".format(maxPhredStrandBias))
        } else {
          None
        }
      })
  }

  /**
   * @param genotype call to evaluate.
   * @param maxDepth The maximum depth to emit variant calls for.
   * @return If this genotype shows excessive depth, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterMaxDepth(genotype: Genotype,
                                       maxDepth: Int): Option[String] = {
    Option(genotype.getReadDepth)
      .flatMap(dp => {
        if (dp > maxDepth) {
          Some("DP>%d".format(maxDepth))
        } else {
          None
        }
      })
  }

  /**
   * @param genotype call to evaluate.
   * @param minDepth The minimum depth to emit variant calls for.
   * @return If this genotype shows insufficient depth, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterMinDepth(genotype: Genotype,
                                       minDepth: Int): Option[String] = {
    Option(genotype.getReadDepth)
      .flatMap(dp => {
        if (dp < minDepth) {
          Some("DP<%d".format(minDepth))
        } else {
          None
        }
      })
  }

  /**
   * @param genotype call to evaluate.
   * @param minRMSMappingQuality The minimum root mean square mapping quality to
   *   allow.
   * @return If this genotype shows low average mapping quality, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterRMSMapQ(genotype: Genotype,
                                      minRMSMappingQuality: Float): Option[String] = {
    Option(genotype.getVariantCallingAnnotations)
      .flatMap(vca => Option(vca.getRmsMapQ))
      .flatMap(rmsMapQ => {
        if (rmsMapQ < minRMSMappingQuality) {
          Some("MQ<%1.1f".format(minRMSMappingQuality))
        } else {
          None
        }
      })
  }

  /**
   * Applies a set of hard filters to a genotype call.
   *
   * @param genotype call to evaluate.
   * @param hardFilters A collection of hard filtering functions to apply.
   * @return Returns the genotype with hard filters applied.
   */
  private[util] def hardFilterGenotype(
    genotype: Genotype,
    hardFilters: Iterable[Genotype => Option[String]]): Genotype = {

    // apply the hard filters
    val failedFilterMsgs = hardFilters.flatMap(filterFn => filterFn(genotype))

    // update the genotype
    updateGenotypeWithFilters(genotype,
      hardFilters.nonEmpty,
      failedFilterMsgs)
  }

  /**
   * @param genotype record to update.
   * @param filtersWereApplied True if we applied filters to this variant.
   * @param failedFilters Collection of messages for hard filters that failed.
   * @return Returns an updated genotype record.
   */
  private[util] def updateGenotypeWithFilters(
    genotype: Genotype,
    filtersWereApplied: Boolean,
    failedFilters: Iterable[String]): Genotype = {

    // do we have variant calling annotations? start building anew
    val vcab = Option(genotype.getVariantCallingAnnotations)
      .fold(VariantCallingAnnotations.newBuilder)(vca => {
        VariantCallingAnnotations.newBuilder(vca)
      })

    // set whether we applied filters, and our filter array
    if (filtersWereApplied) {
      vcab.setFiltersApplied(true)
        .setFiltersPassed(failedFilters.isEmpty)
    }
    if (failedFilters.nonEmpty) {
      vcab.setFiltersFailed(failedFilters.toSeq)
    }

    // replace/add our variant calling annotations, build, and return
    Genotype.newBuilder(genotype)
      .setVariantCallingAnnotations(vcab.build)
      .build
  }

  private def siteIsSnp(genotype: Genotype): Boolean = {
    genotype.getVariant.getReferenceAllele.length == 1 &&
      genotype.getVariant.getAlternateAllele.length == 1
  }

  /**
   * Filters genotype for emission, and applies hard filters.
   *
   * @param genotype call to evaluate.
   * @param minQuality The minimum genotype quality to emit.
   * @param snpFilters Collection of filters to apply to emitted SNPs.
   * @param indelFilters Collection of filters to apply to emitted INDELs.
   * @return If genotype is high enough quality to emit, a hard filtered
   *   genotype.
   */
  private[util] def filterGenotype(
    genotype: Genotype,
    minQuality: Int,
    snpFilters: Iterable[Genotype => Option[String]],
    indelFilters: Iterable[Genotype => Option[String]]): Option[Genotype] = {

    // first, apply emission filters
    val optGenotype = Some(genotype).filter(emitGenotypeFilter(_, minQuality))

    // then, check whether we are a snp or indel and apply hard filters
    optGenotype.map(gt => {
      if (siteIsSnp(gt)) {
        hardFilterGenotype(gt, snpFilters)
      } else {
        hardFilterGenotype(gt, indelFilters)
      }
    })
  }
}
