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
package org.bdgenomics.avocado.cli

import org.apache.spark.SparkContext
import org.bdgenomics.adam.projections.{ AlignmentRecordField, Filter }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.ADAMSaveAnyArgs
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variant.GenotypeRDD
import org.bdgenomics.avocado.genotyping.{
  BiallelicGenotyper => Biallelic,
  DiscoverVariants => Discover,
  TrioCaller
}
import org.bdgenomics.avocado.util.{
  HardFilterGenotypes,
  HardFilterGenotypesArgs,
  PrefilterReads,
  PrefilterReadsArgs
}
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.utils.cli._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

object TrioGenotyper extends BDGCommandCompanion {
  val commandName = "trioGenotyper"
  val commandDescription = "Call variants in a trio under a biallelic model"

  def apply(cmdLine: Array[String]) = {
    new TrioGenotyper(Args4j[TrioGenotyperArgs](cmdLine))
  }
}

class TrioGenotyperArgs extends Args4jBase with ADAMSaveAnyArgs with ParquetArgs with PrefilterReadsArgs with HardFilterGenotypesArgs {
  @Argument(required = true,
    metaVar = "FIRST_PARENT",
    usage = "The ADAM, BAM or SAM file to call",
    index = 0)
  var firstParentPath: String = null
  @Argument(required = true,
    metaVar = "SECOND_PARENT",
    usage = "The ADAM, BAM or SAM file to call",
    index = 1)
  var secondParentPath: String = null
  @Argument(required = true,
    metaVar = "CHILD",
    usage = "The ADAM, BAM or SAM file to call",
    index = 2)
  var childPath: String = null
  @Argument(required = true,
    metaVar = "OUTPUT",
    usage = "Location to write the variant output",
    index = 3)
  var outputPath: String = null
  @Args4jOption(required = false,
    name = "-ploidy",
    usage = "Ploidy to call using. Defaults to diploid (2).")
  var ploidy: Int = 2
  @Args4jOption(required = false,
    name = "-variants_to_call",
    usage = "Optional previously discovered variants to call. Will skip variant discovery.")
  var variantsToCall: String = _
  @Args4jOption(required = false,
    name = "-is_not_grc",
    usage = "True if the genome build is not from the GRC, or does not have GRC chr prefixes")
  var isNotGrc: Boolean = false
  @Args4jOption(required = false,
    name = "-autosomal_only",
    usage = "True if we want to discard the sex chromosomes.")
  var autosomalOnly: Boolean = false
  @Args4jOption(required = false,
    name = "-keep_mitochondrial_chromosome",
    usage = "True if we want to call variants on the mitochondrial chromosome")
  var keepMitochondrialChromosome: Boolean = false
  @Args4jOption(required = false,
    name = "-keep_duplicate_reads",
    usage = "True if we want to include reads that were marked as duplicates")
  var keepDuplicates: Boolean = false
  @Args4jOption(required = false,
    name = "-min_mapping_quality",
    usage = "Minimum read mapping quality to keep. Defaults to MapQ = 10.")
  var minMappingQuality: Int = 10
  @Args4jOption(required = false,
    name = "-keep_non_primary",
    usage = "If set, we keep secondary/supplemental alignments.")
  var keepNonPrimary: Boolean = false
  @Args4jOption(required = false,
    name = "-min_phred_to_discover_variant",
    usage = "Minimum quality needed to discover a variant. Defaults to phred 25.")
  var minPhredForDiscovery: Int = 25
  @Args4jOption(required = false,
    name = "-min_observations_to_discover_variant",
    usage = "Minimum number of times a variant must be seen to be discovered. Defaults to phred 5.")
  var minObservationsForDiscovery: Int = 5
  @Args4jOption(required = false,
    name = "-min_genotype_quality",
    usage = "Minimum quality needed to emit a non-ref genotype. Default is Phred 30.")
  var minQuality: Int = 30
  @Args4jOption(required = false,
    name = "-min_het_snp_quality_by_depth",
    usage = "Minimum heterozygous SNP quality/depth for hard filtering. Default is 2.0. Set negative to ignore filter.")
  var minHetSnpQualityByDepth: Float = 2.0f
  @Args4jOption(required = false,
    name = "-min_hom_snp_quality_by_depth",
    usage = "Minimum homozygous SNP quality/depth for hard filtering. Default is 1.0. Set negative to ignore filter.")
  var minHomSnpQualityByDepth: Float = 1.0f
  @Args4jOption(required = false,
    name = "-min_het_indel_quality_by_depth",
    usage = "Minimum heterozygous INDEL quality/depth for hard filtering. Default is 2.0. Set negative to ignore filter.")
  var minHetIndelQualityByDepth: Float = 2.0f
  @Args4jOption(required = false,
    name = "-min_hom_indel_quality_by_depth",
    usage = "Minimum homozygous INDEL quality/depth for hard filtering. Default is 1.0. Set negative to ignore filter.")
  var minHomIndelQualityByDepth: Float = 1.0f
  @Args4jOption(required = false,
    name = "-min_snp_rms_mapping_quality",
    usage = "Minimum SNP root mean square mapping quality for hard filtering. Default is 30.0. Set negative to ignore filter.")
  var minSnpRMSMappingQuality: Float = 30.0f
  @Args4jOption(required = false,
    name = "-min_indel_rms_mapping_quality",
    usage = "Minimum INDEL root mean square mapping quality for hard filtering. Default is -1.0 (ignored). Set negative to ignore filter.")
  var minIndelRMSMappingQuality: Float = -1.0f
  @Args4jOption(required = false,
    name = "-min_snp_depth",
    usage = "Minimum SNP depth for hard filtering. Default is 10. Set to a negative value to omit.")
  var minSnpDepth: Int = 10
  @Args4jOption(required = false,
    name = "-max_snp_depth",
    usage = "Maximum SNP depth for hard filtering. Default is 300. Set to a negative value to omit.")
  var maxSnpDepth: Int = 200
  @Args4jOption(required = false,
    name = "-min_indel_depth",
    usage = "Minimum INDEL depth for hard filtering. Default is 10. Set to a negative value to omit.")
  var minIndelDepth: Int = 10
  @Args4jOption(required = false,
    name = "-max_indel_depth",
    usage = "Maximum INDEL depth for hard filtering. Default is 300. Set to a negative value to omit.")
  var maxIndelDepth: Int = 200
  @Args4jOption(required = false,
    name = "-min_het_snp_allelic_fraction",
    usage = "Minimum (alt) allelic fraction for calling a het SNP. Default is 0.333. Set to a negative value to omit.")
  var minHetSnpAltAllelicFraction: Float = 0.333f
  @Args4jOption(required = false,
    name = "-max_het_snp_allelic_fraction",
    usage = "Maximum (alt) allelic fraction for calling a het SNP. Default is 0.666. Set to a negative value to omit.")
  var maxHetSnpAltAllelicFraction: Float = 0.666f
  @Args4jOption(required = false,
    name = "-min_hom_snp_allelic_fraction",
    usage = "Minimum (alt) allelic fraction for calling a hom SNP. Default is 0.666. Set to a negative value to omit.")
  var minHomSnpAltAllelicFraction: Float = 0.666f
  @Args4jOption(required = false,
    name = "-min_het_indel_allelic_fraction",
    usage = "Minimum (alt) allelic fraction for calling a het INDEL. Default is 0.333. Set to a negative value to omit.")
  var minHetIndelAltAllelicFraction: Float = 0.333f
  @Args4jOption(required = false,
    name = "-max_het_indel_allelic_fraction",
    usage = "Maximum (alt) allelic fraction for calling a het INDEL. Default is 0.666. Set to a negative value to omit.")
  var maxHetIndelAltAllelicFraction: Float = 0.666f
  @Args4jOption(required = false,
    name = "-min_hom_indel_allelic_fraction",
    usage = "Minimum (alt) allelic fraction for calling a hom INDEL. Default is 0.666. Set to a negative value to omit.")
  var minHomIndelAltAllelicFraction: Float = 0.666f

  // required by HardFilterGenotypesArgs
  var maxSnpPhredStrandBias: Float = -1.0f
  var maxIndelPhredStrandBias: Float = -1.0f

  // required by ADAMSaveAnyArgs
  var sortFastqOutput: Boolean = false
  var asSingleFile: Boolean = false
  var deferMerging: Boolean = false
  var disableFastConcat: Boolean = false
}

class TrioGenotyper(
    protected val args: TrioGenotyperArgs) extends BDGSparkCommand[TrioGenotyperArgs] {

  val companion = TrioGenotyper

  def run(sc: SparkContext) {

    // load reads
    val projection = Some(Filter(AlignmentRecordField.attributes,
      AlignmentRecordField.origQual,
      AlignmentRecordField.recordGroupName))
    val firstParentReads = sc.loadAlignments(args.firstParentPath,
      optProjection = projection)
    val secondParentReads = sc.loadAlignments(args.secondParentPath,
      optProjection = projection)
    val childReads = sc.loadAlignments(args.childPath,
      optProjection = projection)

    // get sample IDs
    val firstParentId = TrioCaller.extractSampleId(firstParentReads)
    val secondParentId = TrioCaller.extractSampleId(secondParentReads)
    val childId = TrioCaller.extractSampleId(childReads)

    // merge reads
    val reads = firstParentReads.union(secondParentReads, childReads)

    // filter reads
    val filteredReads = PrefilterReads(reads, args)

    // were we provided variants? if so, load them and call.
    // else, discover variants and call
    val variants = Option(args.variantsToCall).fold({
      Discover(filteredReads,
        optPhredThreshold = Some(args.minPhredForDiscovery),
        optMinObservations = Some(args.minObservationsForDiscovery))
    })(vPath => {

      // load variants
      sc.loadVariants(vPath)
    })

    val firstParentGenotypes = Biallelic.call(
      PrefilterReads(firstParentReads, args),
      variants,
      args.ploidy)

    val secondParentGenotypes = Biallelic.call(
      PrefilterReads(secondParentReads, args),
      variants,
      args.ploidy)

    val childGenotypes = Biallelic.call(
      PrefilterReads(childReads, args),
      variants,
      args.ploidy)

    /*val genotypes = GenotypeRDD(sc.union(firstParentGenotypes.rdd,
      secondParentGenotypes.rdd,
      childGenotypes.rdd),
      variants.sequences,
      (firstParentGenotypes.samples ++
        secondParentGenotypes.samples ++
        childGenotypes.samples))

    // hard filter the genotypes
    val filteredGenotypes = HardFilterGenotypes(genotypes,
      args,
      filterRefGenotypes = false)

    // trio call the site
    val trioGenotypes = TrioCaller(genotypes,
      firstParentId,
      secondParentId,
      childId)

    // save the variant calls
    trioGenotypes.saveAsParquet(args)*/
  }
}
