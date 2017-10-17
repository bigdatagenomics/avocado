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
import org.bdgenomics.avocado.genotyping.{ BiallelicGenotyper => Biallelic }
import org.bdgenomics.avocado.models.CopyNumberMap
import org.bdgenomics.avocado.util.{
  HardFilterGenotypes,
  HardFilterGenotypesArgs,
  PrefilterReads,
  PrefilterReadsArgs,
  RewriteHets,
  RewriteHetsArgs
}
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.utils.cli._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

object BiallelicGenotyper extends BDGCommandCompanion {
  val commandName = "biallelicGenotyper"
  val commandDescription = "Call variants under a biallelic model"

  def apply(cmdLine: Array[String]) = {
    new BiallelicGenotyper(Args4j[BiallelicGenotyperArgs](cmdLine))
  }
}

class BiallelicGenotyperArgs extends Args4jBase with ADAMSaveAnyArgs with ParquetArgs with PrefilterReadsArgs with HardFilterGenotypesArgs with RewriteHetsArgs {
  @Argument(required = true,
    metaVar = "INPUT",
    usage = "The ADAM, BAM or SAM file to call",
    index = 0)
  var inputPath: String = null
  @Argument(required = true,
    metaVar = "OUTPUT",
    usage = "Location to write the variant output",
    index = 1)
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
    name = "-desired_partition_count",
    usage = "The desired number of partitions to have after the shuffle. Negative values are ignored. Defaults to -1.")
  var desiredPartitionCount: Int = -1
  @Args4jOption(required = false,
    name = "-desired_partition_size",
    usage = "The desired number of reads per partition. Negative values are ignored.")
  var desiredPartitionSize: Int = -1
  @Args4jOption(required = false,
    name = "-desired_max_coverage",
    usage = "The maximum desired number of reads per site. Negative values are ignored. Default is 1000.")
  var desiredMaxCoverage: Int = 1000
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
  @Args4jOption(required = false,
    name = "-disable_het_snp_rewriting",
    usage = "If true, disables rewriting of high allelic fraction het SNPs as hom alt SNPs.")
  var disableHetSnpRewriting: Boolean = false
  @Args4jOption(required = false,
    name = "-disable_het_indel_rewriting",
    usage = "If true, disables rewriting of high allelic fraction het INDELs as hom alt INDELs.")
  var disableHetIndelRewriting: Boolean = false
  @Args4jOption(required = false,
    name = "-score_all_sites",
    usage = "If provided, scores all sites, even non-variant sites. Emits a gVCF styled output.")
  var scoreAllSites = false
  @Args4jOption(required = false,
    name = "-cnvs",
    usage = "Copy number variant calls for this sample.")
  var cnvCalls: String = null
  @Args4jOption(required = false,
    name = "-emit_all_genotypes",
    usage = "If true, emits all genotyped sites. Use if joint calling.")
  var emitAllGenotypes = false

  // required by HardFilterGenotypesArgs
  var maxSnpPhredStrandBias: Float = -1.0f
  var maxIndelPhredStrandBias: Float = -1.0f

  // required by ADAMSaveAnyArgs
  var sortFastqOutput: Boolean = false
  var asSingleFile: Boolean = false
  var deferMerging: Boolean = false
  var disableFastConcat: Boolean = false
}

class BiallelicGenotyper(
    protected val args: BiallelicGenotyperArgs) extends BDGSparkCommand[BiallelicGenotyperArgs] {

  val companion = BiallelicGenotyper

  def run(sc: SparkContext) {

    // load reads
    val projection = Some(Filter(AlignmentRecordField.attributes,
      AlignmentRecordField.origQual,
      AlignmentRecordField.recordGroupName))
    val reads = sc.loadAlignments(args.inputPath,
      optProjection = projection)
    val samples = reads.recordGroups.recordGroups.map(_.sample).toSet
    require(samples.size <= 1,
      "Saw more than one sample (%s) attached to input.".format(samples.mkString(", ")))

    // filter reads
    val filteredReads = PrefilterReads(reads, args)

    // has the user optionally requested a number of partitions?
    val optDesiredPartitionCount = Option(args.desiredPartitionCount)
      .filter(_ >= 1)

    // has the user optionally requested a partition size?
    val optDesiredPartitionSize = Option(args.desiredPartitionSize)
      .filter(_ >= 1)

    // has the user asked us to filter the reads?
    val optDesiredMaxCoverage = Option(args.desiredMaxCoverage)
      .filter(_ >= 1)

    // has the user provided copy number variant calls?
    val copyNumber = Option(args.cnvCalls)
      .fold(CopyNumberMap.empty(args.ploidy))(p => {
        // select cnv calls where the source is a known sample ID
        val features = sc.loadFeatures(p)
          .transform(_.filter(f => samples(f.getSource)))
        CopyNumberMap(args.ploidy, features)
      })

    // were we provided variants? if so, load them and call.
    // else, discover variants and call
    val genotypes = Option(args.variantsToCall).fold({
      Biallelic.discoverAndCall(filteredReads,
        copyNumber,
        args.scoreAllSites,
        optDesiredPartitionCount = optDesiredPartitionCount,
        optPhredThreshold = Some(args.minPhredForDiscovery),
        optMinObservations = Some(args.minObservationsForDiscovery),
        optDesiredPartitionSize = optDesiredPartitionSize,
        optDesiredMaxCoverage = optDesiredMaxCoverage)
    })(vPath => {

      // load variants
      val variants = sc.loadVariants(vPath)

      Biallelic.call(filteredReads,
        variants,
        copyNumber,
        args.scoreAllSites,
        optDesiredPartitionCount = optDesiredPartitionCount,
        optDesiredPartitionSize = optDesiredPartitionSize,
        optDesiredMaxCoverage = optDesiredMaxCoverage)
    })

    // hard filter the genotypes
    val filteredGenotypes = HardFilterGenotypes(RewriteHets(genotypes, args),
      args,
      filterRefGenotypes = !args.scoreAllSites,
      emitAllGenotypes = args.emitAllGenotypes)

    // save the variant calls
    filteredGenotypes.saveAsParquet(args)
  }
}
