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

import htsjdk.variant.vcf.{ VCFFilterHeaderLine, VCFHeaderLine }
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

  /**
   * The minimum allele fraction for the alternate allele in a het SNP call.
   *
   * Set to a negative value to omit.
   */
  var minHetSnpAltAllelicFraction: Float

  /**
   * The maximum allele fraction for the alternate allele in a het SNP call.
   *
   * Set to a negative value to omit.
   */
  var maxHetSnpAltAllelicFraction: Float

  /**
   * The minimum allele fraction for the alternate allele in a hom SNP call.
   *
   * Set to a negative value to omit.
   */
  var minHomSnpAltAllelicFraction: Float

  /**
   * The minimum allele fraction for the alternate allele in a het INDEL call.
   *
   * Set to a negative value to omit.
   */
  var minHetIndelAltAllelicFraction: Float

  /**
   * The maximum allele fraction for the alternate allele in a het INDEL call.
   *
   * Set to a negative value to omit.
   */
  var maxHetIndelAltAllelicFraction: Float

  /**
   * The minimum allele fraction for the alternate allele in a hom INDEL call.
   *
   * Set to a negative value to omit.
   */
  var minHomIndelAltAllelicFraction: Float
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

    // add headers
    val filterHeaders: Seq[VCFHeaderLine] = Seq(
      Option(args.minHetSnpQualityByDepth).filter(_ > 0.0)
        .map(qd => new VCFFilterHeaderLine("HETSNPQD",
          "Quality by depth was below %f for a heterozygous SNP.".format(qd))),
      Option(args.minHomSnpQualityByDepth).filter(_ > 0.0)
        .map(qd => new VCFFilterHeaderLine("HOMSNPQD",
          "Quality by depth was below %f for a homozygous SNP.".format(qd))),
      Option(args.minHetIndelQualityByDepth).filter(_ > 0.0)
        .map(qd => new VCFFilterHeaderLine("HETINDELQD",
          "Quality by depth was below %f for a heterozygous INDEL.".format(qd))),
      Option(args.minHomIndelQualityByDepth).filter(_ > 0.0)
        .map(qd => new VCFFilterHeaderLine("HOMINDELQD",
          "Quality by depth was below %f for a homozygous INDEL.".format(qd))),
      Option(args.maxSnpPhredStrandBias).filter(_ > 0.0)
        .map(fs => new VCFFilterHeaderLine("SNPFS",
          "Phred Fisher scored strand bias was above %f for a SNP.".format(fs))),
      Option(args.maxIndelPhredStrandBias).filter(_ > 0.0)
        .map(fs => new VCFFilterHeaderLine("INDELFS",
          "Phred Fisher scored strand bias was above %f for a INDEL.".format(fs))),
      Option(args.minSnpRMSMappingQuality).filter(_ > 0.0)
        .map(mq => new VCFFilterHeaderLine("SNPMQ",
          "RMS mapping quality was below %f for a SNP.".format(mq))),
      Option(args.minIndelRMSMappingQuality).filter(_ > 0.0)
        .map(mq => new VCFFilterHeaderLine("INDELMQ",
          "RMS mapping quality was below %f for a INDEL.".format(mq))),
      Option(args.minSnpDepth).filter(_ > 0)
        .map(dp => new VCFFilterHeaderLine("SNPMINDP",
          "Read depth was below %d for a SNP.".format(dp))),
      Option(args.minIndelDepth).filter(_ > 0)
        .map(dp => new VCFFilterHeaderLine("INDELMINDP",
          "Read depth was below %d for a INDEL.".format(dp))),
      Option(args.maxSnpDepth).filter(_ > 0)
        .map(dp => new VCFFilterHeaderLine("SNPMAXDP",
          "Read depth was above %d for a SNP.".format(dp))),
      Option(args.maxIndelDepth).filter(_ > 0)
        .map(dp => new VCFFilterHeaderLine("INDELMAXDP",
          "Read depth was above %d for a INDEL.".format(dp))),
      Option(args.minHetSnpAltAllelicFraction).filter(_ > 0)
        .map(af => new VCFFilterHeaderLine("HETSNPMINAF",
          "Allelic fraction was below %f for a het SNP.".format(af))),
      Option(args.maxHetSnpAltAllelicFraction).filter(_ > 0)
        .map(af => new VCFFilterHeaderLine("HETSNPMAXAF",
          "Allelic fraction was above %f for a het SNP.".format(af))),
      Option(args.minHomSnpAltAllelicFraction).filter(_ > 0)
        .map(af => new VCFFilterHeaderLine("HOMSNPMINAF",
          "Allelic fraction was below %f for a hom SNP.".format(af))),
      Option(args.minHetSnpAltAllelicFraction).filter(_ > 0)
        .map(af => new VCFFilterHeaderLine("HETINDELMINAF",
          "Allelic fraction was below %f for a het INDEL.".format(af))),
      Option(args.maxHetSnpAltAllelicFraction).filter(_ > 0)
        .map(af => new VCFFilterHeaderLine("HETINDELMAXAF",
          "Allelic fraction was above %f for a het INDEL.".format(af))),
      Option(args.minHomSnpAltAllelicFraction).filter(_ > 0)
        .map(af => new VCFFilterHeaderLine("HOMINDELMINAF",
          "Allelic fraction was below %f for a hom INDEL.".format(af))))
      .flatten

    // flat map the filters over the genotype rdd
    val minQuality = args.minQuality
    grdd.transform(rdd => {
      rdd.flatMap(filterGenotype(_,
        minQuality,
        snpFilters,
        indelFilters))
    }).copy(headerLines = (grdd.headerLines ++ filterHeaders).distinct)
  }

  /**
   * @param args The hard filter configuration to apply.
   * @return A collection of hard filtering functions to apply to SNPs.
   */
  private[util] def buildSnpHardFilters(
    args: HardFilterGenotypesArgs): Iterable[Genotype => Option[String]] = {
    val optHetQdFilter: Option[Genotype => Option[String]] = Option(args.minHetSnpQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, "HETSNPQD", hom = false))
    val optHomQdFilter: Option[Genotype => Option[String]] = Option(args.minHomSnpQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, "HOMSNPQD", hom = true))
    val optMqFilter: Option[Genotype => Option[String]] = Option(args.minSnpRMSMappingQuality)
      .filter(_ > 0.0)
      .map(minMq => hardFilterRMSMapQ(_, minMq, "SNPMQ"))
    val optFsFilter: Option[Genotype => Option[String]] = Option(args.maxSnpPhredStrandBias)
      .filter(_ > 0.0)
      .map(minFs => hardFilterStrandBias(_, minFs, "SNPFS"))
    val optMinDpFilter: Option[Genotype => Option[String]] = Option(args.minSnpDepth)
      .filter(_ > 0)
      .map(minDp => hardFilterMinDepth(_, minDp, "SNPMINDP"))
    val optMaxDpFilter: Option[Genotype => Option[String]] = Option(args.maxSnpDepth)
      .filter(_ > 0)
      .map(maxDp => hardFilterMaxDepth(_, maxDp, "SNPMAXDP"))
    val optMinAfHetFilter: Option[Genotype => Option[String]] = Option(args.minHetSnpAltAllelicFraction)
      .filter(_ > 0)
      .map(minAf => hardFilterMinAllelicFraction(_, minAf, "HETSNPMINAF", hom = false))
    val optMaxAfHetFilter: Option[Genotype => Option[String]] = Option(args.maxHetSnpAltAllelicFraction)
      .filter(_ > 0)
      .map(maxAf => hardFilterMinAllelicFraction(_, maxAf, "HETSNPMAXAF", hom = false))
    val optMinAfHomFilter: Option[Genotype => Option[String]] = Option(args.minHomSnpAltAllelicFraction)
      .filter(_ > 0)
      .map(minAf => hardFilterMinAllelicFraction(_, minAf, "HOMSNPMINAF", hom = true))

    Iterable(optHetQdFilter,
      optHomQdFilter,
      optMqFilter,
      optFsFilter,
      optMinDpFilter,
      optMaxDpFilter,
      optMinAfHetFilter,
      optMaxAfHetFilter,
      optMinAfHomFilter).flatten
  }

  /**
   * @param args The hard filter configuration to apply.
   * @return A collection of hard filtering functions to apply to INDELs.
   */
  private[util] def buildIndelHardFilters(
    args: HardFilterGenotypesArgs): Iterable[Genotype => Option[String]] = {
    val optHetQdFilter: Option[Genotype => Option[String]] = Option(args.minHetIndelQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, "HETINDELQD", hom = false))
    val optHomQdFilter: Option[Genotype => Option[String]] = Option(args.minHomIndelQualityByDepth)
      .filter(_ > 0.0)
      .map(minQd => hardFilterQualityByDepth(_, minQd, "HOMINDELQD", hom = true))
    val optMqFilter: Option[Genotype => Option[String]] = Option(args.minIndelRMSMappingQuality)
      .filter(_ > 0.0)
      .map(minMq => hardFilterRMSMapQ(_, minMq, "INDELMQ"))
    val optFsFilter: Option[Genotype => Option[String]] = Option(args.maxIndelPhredStrandBias)
      .filter(_ > 0.0)
      .map(minFs => hardFilterStrandBias(_, minFs, "INDELFS"))
    val optMinDpFilter: Option[Genotype => Option[String]] = Option(args.minIndelDepth)
      .filter(_ > 0)
      .map(minDp => hardFilterMinDepth(_, minDp, "INDELMINDP"))
    val optMaxDpFilter: Option[Genotype => Option[String]] = Option(args.maxIndelDepth)
      .filter(_ > 0)
      .map(maxDp => hardFilterMaxDepth(_, maxDp, "INDELMAXDP"))
    val optMinAfHetFilter: Option[Genotype => Option[String]] = Option(args.minHetIndelAltAllelicFraction)
      .filter(_ > 0)
      .map(minAf => hardFilterMinAllelicFraction(_, minAf, "HETINDELMINAF", hom = false))
    val optMaxAfHetFilter: Option[Genotype => Option[String]] = Option(args.maxHetIndelAltAllelicFraction)
      .filter(_ > 0)
      .map(maxAf => hardFilterMinAllelicFraction(_, maxAf, "HETINDELMAXAF", hom = false))
    val optMinAfHomFilter: Option[Genotype => Option[String]] = Option(args.minHomIndelAltAllelicFraction)
      .filter(_ > 0)
      .map(minAf => hardFilterMinAllelicFraction(_, minAf, "HOMINDELMINAF", hom = true))

    Iterable(optHetQdFilter,
      optHomQdFilter,
      optMqFilter,
      optFsFilter,
      optMinDpFilter,
      optMaxDpFilter,
      optMinAfHetFilter,
      optMaxAfHetFilter,
      optMinAfHomFilter).flatten
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
   * @param msg Filter message to use.
   * @param hom If true, this filter should be applied to homozygous variants.
   * @return If this genotype has a low quality for the depth of the site, we
   *   return an optional filter string. Else, we return a None.
   */
  private[util] def hardFilterQualityByDepth(genotype: Genotype,
                                             minQualityByDepth: Float,
                                             msg: String,
                                             hom: Boolean = false): Option[String] = {
    val gtIsHom = genotype.getAlleles.forall(_ == GenotypeAllele.ALT)
    if ((!gtIsHom && !hom) || (gtIsHom && hom)) {
      Option(genotype.getReadDepth)
        .flatMap(depth => {
          Option(genotype.getGenotypeQuality)
            .map(gq => gq.toFloat / depth.toFloat)
        }).flatMap(qualByDepth => {
          if (qualByDepth < minQualityByDepth) {
            Some(msg)
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
   * @param msg Filter message to use.
   * @return If this genotype shows significant strand bias, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterStrandBias(genotype: Genotype,
                                         maxPhredStrandBias: Float,
                                         msg: String): Option[String] = {
    Option(genotype.getVariantCallingAnnotations)
      .flatMap(vca => Option(vca.getFisherStrandBiasPValue))
      .flatMap(fs => {
        if (fs > maxPhredStrandBias) {
          Some(msg)
        } else {
          None
        }
      })
  }

  /**
   * @param genotype call to evaluate.
   * @param maxDepth The maximum depth to emit variant calls for.
   * @param msg Filter message to use.
   * @return If this genotype shows excessive depth, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterMaxDepth(genotype: Genotype,
                                       maxDepth: Int,
                                       msg: String): Option[String] = {
    Option(genotype.getReadDepth)
      .flatMap(dp => {
        if (dp > maxDepth) {
          Some(msg)
        } else {
          None
        }
      })
  }

  /**
   * @param genotype call to evaluate.
   * @param minDepth The minimum depth to emit variant calls for.
   * @param msg Filter message to use.
   * @return If this genotype shows insufficient depth, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterMinDepth(genotype: Genotype,
                                       minDepth: Int,
                                       msg: String): Option[String] = {
    Option(genotype.getReadDepth)
      .flatMap(dp => {
        if (dp < minDepth) {
          Some(msg)
        } else {
          None
        }
      })
  }

  /**
   * @param genotype call to evaluate.
   * @param minRMSMappingQuality The minimum root mean square mapping quality to
   *   allow.
   * @param msg Filter message to use.
   * @return If this genotype shows low average mapping quality, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterRMSMapQ(genotype: Genotype,
                                      minRMSMappingQuality: Float,
                                      msg: String): Option[String] = {
    Option(genotype.getVariantCallingAnnotations)
      .flatMap(vca => Option(vca.getRmsMapQ))
      .flatMap(rmsMapQ => {
        if (rmsMapQ < minRMSMappingQuality) {
          Some(msg)
        } else {
          None
        }
      })
  }

  /**
   * @param genotype call to evaluate.
   * @return Returns the alleleic fraction for the alt allele.
   */
  private[util] def optAlleleFraction(genotype: Genotype): Option[Float] = {
    (Option(genotype.getReadDepth), Option(genotype.getAlternateReadDepth)) match {
      case (Some(dp), Some(alt)) => Some(alt.toFloat / dp.toFloat)
      case _                     => None
    }
  }

  /**
   * @param genotype call to evaluate.
   * @param minAlleleFraction The minimum allele fraction to allow.
   * @param msg Filter message to use.
   * @param hom If true, only consider hom sites; else only consider hets.
   * @return If this genotype shows low allele fraction, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterMinAllelicFraction(genotype: Genotype,
                                                 minAlleleFraction: Float,
                                                 msg: String,
                                                 hom: Boolean = false): Option[String] = {
    val gtIsHom = genotype.getAlleles.forall(_ == GenotypeAllele.ALT)
    if ((!gtIsHom && !hom) || (gtIsHom && hom)) {
      optAlleleFraction(genotype)
        .flatMap(af => {
          if (af <= minAlleleFraction) {
            Some(msg)
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
   * @param maxAlleleFraction The maximum allele fraction to allow.
   * @param msg Filter message to use.
   * @return If this genotype shows low allele fraction, we return an
   *   optional filter string. Else, we return a None.
   */
  private[util] def hardFilterMaxAllelicFraction(genotype: Genotype,
                                                 maxAlleleFraction: Float,
                                                 msg: String): Option[String] = {
    val gtIsHom = genotype.getAlleles.forall(_ == GenotypeAllele.ALT)
    if (gtIsHom) {
      None
    } else {
      optAlleleFraction(genotype)
        .flatMap(af => {
          if (af > maxAlleleFraction) {
            Some(msg)
          } else {
            None
          }
        })
    }
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
