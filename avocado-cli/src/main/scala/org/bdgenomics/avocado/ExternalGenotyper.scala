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
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.{ ADAMSaveAnyArgs, InFormatterCompanion }
import org.bdgenomics.adam.rdd.read.{
  AlignmentRecordRDD,
  BAMInFormatter,
  SAMInFormatter
}
import org.bdgenomics.adam.rdd.variation.{
  VariantContextRDD,
  VCFOutFormatter
}
import org.bdgenomics.avocado.util.{ PrefilterReads, PrefilterReadsArgs }
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.utils.cli._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

object ExternalGenotyper extends BDGCommandCompanion {
  val commandName = "externalGenotyper"
  val commandDescription = "Call variants by piping to an external tool"

  def apply(cmdLine: Array[String]) = {
    new ExternalGenotyper(Args4j[ExternalGenotyperArgs](cmdLine))
  }
}

class ExternalGenotyperArgs extends Args4jBase with ADAMSaveAnyArgs with ParquetArgs with PrefilterReadsArgs {
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
  @Argument(required = true,
    metaVar = "COMMAND",
    usage = "Piped command to run",
    index = 2)
  var command: String = null
  @Args4jOption(required = false,
    name = "-single",
    usage = "Saves OUTPUT as single file")
  var asSingleFile: Boolean = false
  @Args4jOption(required = false,
    name = "-defer_merging",
    usage = "Defers merging single file output")
  var deferMerging: Boolean = false
  @Args4jOption(required = false,
    name = "-files_to_provide",
    usage = "Comma separated list of files to make accessible to the piped command")
  var filesToProvide: String = _
  @Args4jOption(required = false,
    name = "-env_to_set",
    usage = "Comma separated ENVAR=VALUE pairs to set in each piped command")
  var envToSet: String = _
  @Args4jOption(required = false,
    name = "-save_as_vcf",
    usage = "Save the output as VCF.")
  var saveAsVcf: Boolean = false
  @Args4jOption(required = false,
    name = "-pipe_sam",
    usage = "Pipes the input as SAM. If not set, pipes BAM instead.")
  var pipeSam: Boolean = false
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

  // required by ADAMSaveAnyArgs
  var sortFastqOutput: Boolean = false
}

class ExternalGenotyper(
    protected val args: ExternalGenotyperArgs) extends BDGSparkCommand[ExternalGenotyperArgs] {

  val companion = ExternalGenotyper

  def run(sc: SparkContext) {

    // load reads
    val reads = sc.loadAlignments(args.inputPath)

    // filter reads
    val filteredReads = PrefilterReads(reads, args)

    // split files
    val fileList = Option(args.filesToProvide).getOrElse("")
      .split(",")
      .filter(_.length > 0)

    // split environment variables
    val envMap = Option(args.envToSet).getOrElse("")
      .split(",")
      .filter(_.length > 0)
      .map(s => {
        val (envar, value) = s.span(_ != '=')
        (envar -> value.tail)
      }).toMap

    // append the file $ to the command, if any
    val cmd = if (fileList.nonEmpty) {
      "%s %s".format(args.command, fileList.indices
        .map(s => "$%d".format(s))
        .mkString(" "))
    } else {
      args.command
    }

    // set up pipe formatters and call variants through a pipe
    implicit val uFormatter = new VCFOutFormatter
    val variants: VariantContextRDD = if (args.pipeSam) {
      implicit val tFormatter = SAMInFormatter
      filteredReads.pipe(cmd,
        files = fileList,
        environment = envMap)
    } else {
      implicit val tFormatter = BAMInFormatter
      filteredReads.pipe(cmd,
        files = fileList,
        environment = envMap)
    }

    // save the variant calls
    if (args.saveAsVcf) {
      variants.saveAsVcf(args.outputPath,
        asSingleFile = args.asSingleFile,
        sortOnSave = false)
    } else {
      variants.toGenotypeRDD
        .saveAsParquet(args)
    }
  }
}
