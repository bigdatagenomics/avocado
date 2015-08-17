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

import java.nio.file.Files
import org.apache.commons.configuration.HierarchicalConfiguration
import org.apache.commons.configuration.plist.PropertyListConfiguration
import org.apache.hadoop.mapreduce.Job
import org.apache.spark.rdd.RDD
import org.apache.spark.{ SparkContext, Logging }
import org.kohsuke.args4j.{ Option => option, Argument }
import org.bdgenomics.adam.models.{ VariantContext, ReferenceRegion }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.discovery.Explore
import org.bdgenomics.avocado.genotyping.CallGenotypes
import org.bdgenomics.avocado.input.Input
import org.bdgenomics.avocado.models.Observation
import org.bdgenomics.avocado.preprocessing.Preprocessor
import org.bdgenomics.avocado.postprocessing.Postprocessor
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.{
  Variant,
  AlignmentRecord,
  NucleotideContigFragment,
  Genotype
}
import org.bdgenomics.utils.cli.{
  BDGSparkCommand,
  BDGCommandCompanion,
  ParquetArgs,
  Args4j,
  Args4jBase
}
import org.bdgenomics.utils.instrumentation._

object Avocado extends BDGCommandCompanion {

  val commandName = "Avocado"
  val commandDescription = "Call variants using avocado and the ADAM preprocessing pipeline."

  def apply(args: Array[String]) = {
    new Avocado(Args4j[AvocadoArgs](args))
  }
}

class AvocadoArgs extends Args4jBase with ParquetArgs {
  @Argument(metaVar = "READS", required = true, usage = "ADAM read-oriented data", index = 0)
  var readInput: String = _

  @Argument(metaVar = "REFERENCE", required = true, usage = "ADAM or FASTA reference genome data", index = 1)
  var referenceInput: String = _

  @Argument(metaVar = "VARIANTS", required = true, usage = "ADAM variant output", index = 2)
  var variantOutput: String = _

  @Argument(metaVar = "CONFIG", required = true, usage = "avocado configuration file", index = 3)
  var configFile: String = _

  @option(name = "-debug", usage = "If set, prints a higher level of debug output.")
  var debug = false

  @option(required = false, name = "-fragment_length", usage = "Sets maximum fragment length. Default value is 10,000. Values greater than 1e9 should be avoided.")
  var fragmentLength: Long = 10000L
}

class Avocado(protected val args: AvocadoArgs) extends BDGSparkCommand[AvocadoArgs] with Logging {

  // companion object to this class - needed for BDGCommand framework
  val companion = Avocado

  // get config off classpath and load into a temp file...
  val stream = Thread.currentThread.getContextClassLoader.getResourceAsStream(args.configFile)
  val tempPath = Files.createTempDirectory("config")
  val tempFilePath = tempPath.resolve("temp.properties")
  Files.copy(stream, tempFilePath)

  // load config
  val config: HierarchicalConfiguration = new PropertyListConfiguration(tempFilePath.toFile)

  val preprocessorNames = getStringArrayFromConfig("preprocessorNames")
  val preprocessorAlgorithms = getStringArrayFromConfig("preprocessorAlgorithms")
  assert(preprocessorNames.length == preprocessorAlgorithms.length,
    "Must have a 1-to-1 mapping between preprocessor names and algorithms.")
  val preprocessingStagesZippedWithNames = preprocessorNames.zip(preprocessorAlgorithms)

  val explorerName = config.getString("explorerName")
  val explorerAlgorithm = config.getString("explorerAlgorithm")

  val genotyperName = config.getString("genotyperName")
  val genotyperAlgorithm = config.getString("genotyperAlgorithm")

  val postprocessorNames = getStringArrayFromConfig("postprocessorNames")
  val postprocessorAlgorithms = getStringArrayFromConfig("postprocessorAlgorithms")
  assert(postprocessorNames.length == postprocessorAlgorithms.length,
    "Must have a 1-to-1 mapping between postprocessor names and algoritms.")
  val postprocessorsZipped = postprocessorNames.zip(postprocessorAlgorithms)

  val debug = args.debug

  private def getStringArrayFromConfig(name: String): Array[String] = {
    config.getStringArray(name).map(_.toString)
  }

  /**
   * Applies several pre-processing steps to the read pipeline. Currently, these are the default
   * steps in the ADAM processing pipeline.
   *
   * @param reads RDD of reads to process.
   * @return RDD containing reads that have been sorted and deduped.
   */
  def preProcessReads(reads: RDD[AlignmentRecord]): RDD[AlignmentRecord] = PreprocessReads.time {
    var processedReads = reads //.cache

    if (debug) {
      log.info("avocado: Preprocessing " + processedReads.count + " reads.")
    }

    // loop over preprocessing stages and apply
    preprocessingStagesZippedWithNames.foreach(p => {
      val (stageName, stageAlgorithm) = p

      log.info("avocado: Running " + stageName)

      // run this preprocessing stage
      processedReads = Preprocessor(processedReads, stageName, stageAlgorithm, config)
    })

    // return processed reads
    processedReads
  }

  /**
   * Applies variant calling algorithms to reads and pileups. Reduces down and returns called variants.
   *
   * @param reads
   * @param stats
   * @return Joined output of variant calling algorithms.
   */
  def callVariants(reads: RDD[AlignmentRecord], stats: AvocadoConfigAndStats): RDD[VariantContext] = CallVariants.time {
    val discoveries: RDD[Observation] = Explore(explorerAlgorithm,
      explorerName,
      reads,
      stats,
      config)
    CallGenotypes(genotyperAlgorithm,
      genotyperName,
      discoveries,
      stats,
      config)
  }

  /**
   * Applies variant post-processing methods to called variants. Post-processing can
   * include methods which modify the information in variant calls, or alternatively,
   * methods that filter out spurious variant calls.
   *
   * @param variants RDD of variants to process.
   * @return Post-processed variants.
   */
  def postProcessVariants(variants: RDD[VariantContext], stats: AvocadoConfigAndStats): RDD[VariantContext] = PostprocessVariants.time {
    var rdd = variants

    // loop over post processing steps
    postprocessorsZipped.foreach(p => {
      val (ppStageName, ppAlgorithm) = p

      rdd = Postprocessor(rdd, ppStageName, ppAlgorithm, stats, config)
    })

    rdd
  }

  /**
   * Main method. Implements body of variant caller. SparkContext and Hadoop Job are provided
   * by the ADAMSparkCommand shell.
   *
   * @param sc SparkContext for RDDs.
   * @param job Hadoop Job container for file I/O.
   */
  def run(sc: SparkContext) {

    log.info("Starting avocado...")

    // load in reference from ADAM file
    val reference: RDD[NucleotideContigFragment] = LoadContigs.time {
      sc.loadSequence(args.referenceInput, fragmentLength = args.fragmentLength)
    }

    log.info("Loading reads in from " + args.readInput)
    // load in reads from ADAM file
    val reads: RDD[AlignmentRecord] = LoadReads.time {
      Input(sc, args.readInput, reference, config)
    }

    // create stats/config item
    val stats = new AvocadoConfigAndStats(sc, args.debug, reads, reference)

    // apply read translation steps
    log.info("Processing reads.")
    val cleanedReads = preProcessReads(reads)

    // call variants on filtered reads and pileups
    log.info("Calling variants.")
    val calledVariants = callVariants(cleanedReads, stats)

    // post process variants
    log.info("Post-processing variants.")
    val processedGenotypes: RDD[Genotype] = postProcessVariants(calledVariants, stats).flatMap(variantContext => variantContext.genotypes)

    // save variants to output file
    log.info("Writing calls to disk.")
    SaveVariants.time {
      processedGenotypes.adamParquetSave(args.variantOutput,
        args.blockSize,
        args.pageSize,
        args.compressionCodec,
        args.disableDictionaryEncoding)
    }
  }
}
