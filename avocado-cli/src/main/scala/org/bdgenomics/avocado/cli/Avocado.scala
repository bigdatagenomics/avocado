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

import org.apache.commons.configuration.HierarchicalConfiguration
import org.apache.commons.configuration.plist.PropertyListConfiguration
import org.apache.hadoop.mapreduce.Job
import org.apache.spark.rdd.RDD
import org.apache.spark.{ SparkContext, Logging }
import org.kohsuke.args4j.{ Option => option, Argument }
import org.bdgenomics.formats.avro.{
  ADAMVariant,
  ADAMRecord,
  ADAMNucleotideContigFragment,
  ADAMGenotype
}
import org.bdgenomics.adam.cli.{
  ADAMSparkCommand,
  ADAMCommandCompanion,
  ParquetArgs,
  SparkArgs,
  Args4j,
  Args4jBase
}
import org.bdgenomics.adam.models.{ ADAMVariantContext, ReferenceRegion }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.calls.{ VariantCall, VariantCaller }
import org.bdgenomics.avocado.input.Input
import org.bdgenomics.avocado.partitioners.Partitioner
import org.bdgenomics.avocado.preprocessing.Preprocessor
import org.bdgenomics.avocado.postprocessing.Postprocessor
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

object Avocado extends ADAMCommandCompanion {

  val commandName = "Avocado"
  val commandDescription = "Call variants using avocado and the ADAM preprocessing pipeline."

  def apply(args: Array[String]) = {
    new Avocado(Args4j[AvocadoArgs](args))
  }
}

class AvocadoArgs extends Args4jBase with ParquetArgs with SparkArgs {
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

class Avocado(protected val args: AvocadoArgs) extends ADAMSparkCommand[AvocadoArgs] with Logging {

  // companion object to this class - needed for ADAMCommand framework
  val companion = Avocado

  // get config
  val config: HierarchicalConfiguration = new PropertyListConfiguration(args.configFile)

  val preprocessorNames = getStringArrayFromConfig("preprocessorNames")
  val preprocessorAlgorithms = getStringArrayFromConfig("preprocessorAlgorithms")
  assert(preprocessorNames.length == preprocessorAlgorithms.length,
    "Must have a 1-to-1 mapping between preprocessor names and algorithms.")
  val preprocessingStagesZippedWithNames = preprocessorNames.zip(preprocessorAlgorithms)

  val callNames = getStringArrayFromConfig("callNames")
  val callAlgorithms = getStringArrayFromConfig("callAlgorithms")
  assert(callNames.length == callAlgorithms.length,
    "Must have a 1-to-1 mapping between variant call names and algorithms.")
  val callZipped = callNames.zip(callAlgorithms)

  val partitionerNames = getStringArrayFromConfig("partitionerNames")
  val partitionerAlgorithms = getStringArrayFromConfig("partitionerAlgorithms")
  assert(partitionerNames.length == partitionerAlgorithms.length,
    "Must have a 1-to-1 mapping between partitioner names and algorithms.")
  val partitionsZipped = partitionerNames.zip(partitionerAlgorithms)

  val postprocessorNames = getStringArrayFromConfig("postprocessorNames")
  val postprocessorAlgorithms = getStringArrayFromConfig("postprocessorAlgorithms")
  assert(postprocessorNames.length == postprocessorAlgorithms.length,
    "Must have a 1-to-1 mapping between postprocessor names and algoritms.")
  val postprocessorsZipped = postprocessorNames.zip(postprocessorAlgorithms)

  assert(callZipped.length == partitionsZipped.length,
    "Must have a 1-to-1 mapping between partitioners and calls.")

  val debug = args.debug

  private def getStringArrayFromConfig(name: String): Array[String] = {
    config.getStringArray(name).map(_.toString)
  }

  /**
   * Assigns reads to a variant calling algorithm.
   *
   * @param reads An RDD of raw reads.
   * @return A map that maps RDDs of reads to variant calling algorithms.
   */
  def partitionReads(reads: RDD[ADAMRecord],
                     stats: AvocadoConfigAndStats): Map[VariantCall, RDD[ADAMRecord]] = {

    var rdd = reads.keyBy(r => ReferenceRegion(r))
      .filter(kv => kv._1.isDefined)
      .map(kv => (kv._1.get, kv._2))
    //.cache()

    var callsets = List[(VariantCall, RDD[ADAMRecord])]()
    val partitioners = partitionsZipped.map(p => Partitioner(reads, config, p._1, p._2, stats))
    val partitionersZippedWithCalls = partitioners.zip(callZipped)

    // loop over specified partitioner algorithms
    partitionersZippedWithCalls.foreach(kv => {
      // extract partitioner and call
      val (p, (callName, callAlg)) = kv

      log.info("Partition by: " + p.companion.partitionerName)

      // generate partition set from current read set
      val partitions = p.computePartitions()

      // generate variant caller
      val call = VariantCaller(callName, callAlg, stats, config, partitions)

      // filter reads that are in the partition into a new call set and add the call
      val filteredReads = rdd.filter(kv => partitions.isInSet(kv._1))
        .map(kv => kv._2)
      if (debug) {
        log.info("avocado: have " + filteredReads.count + " reads in " + p.companion.partitionerName)
      }
      callsets = (call, filteredReads) :: callsets

      // filter out reads that are wholly contained in the last partitoner
      rdd = rdd.filter(kv => partitions.isOutsideOfSet(kv._1))
      if (debug) {
        log.info("avocado: have " + rdd.count() + " reads left.")
      }
    })

    // filter out variant calls that are not actually callable
    callsets.toMap.filterKeys(k => k.isCallable)
  }

  /**
   * Applies several pre-processing steps to the read pipeline. Currently, these are the default
   * steps in the ADAM processing pipeline.
   *
   * @param reads RDD of reads to process.
   * @return RDD containing reads that have been sorted and deduped.
   */
  def preProcessReads(reads: RDD[ADAMRecord]): RDD[ADAMRecord] = {
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
   * @param callsets Map containing read calling algorithm and RDD record pairs.
   * @return Joined output of variant calling algorithms.
   */
  def callVariants(callsets: Map[VariantCall, RDD[ADAMRecord]]): RDD[ADAMVariantContext] = {
    // callset map must not be empty
    assert(!callsets.isEmpty)

    // apply calls and take the union of variants called
    callsets.map(pair => {
      // get call and rdd pair
      val (call, rdd) = pair

      if (debug) {
        log.info("avocado: Running " + call.companion.callName + " on " + rdd.count + " reads.")
      }

      // apply call
      call.call(rdd)
    }).reduce(_ ++ _)
  }

  /**
   * Applies variant post-processing methods to called variants. Post-processing can
   * include methods which modify the information in variant calls, or alternatively,
   * methods that filter out spurious variant calls.
   *
   * @param variants RDD of variants to process.
   * @return Post-processed variants.
   */
  def postProcessVariants(variants: RDD[ADAMVariantContext], stats: AvocadoConfigAndStats): RDD[ADAMVariantContext] = {
    var rdd = variants //.cache()

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
  def run(sc: SparkContext, job: Job) {

    log.info("Starting avocado...")

    // load in reference from ADAM file
    val reference: RDD[ADAMNucleotideContigFragment] = sc.adamSequenceLoad(args.referenceInput, args.fragmentLength)

    log.info("Loading reads in from " + args.readInput)
    // load in reads from ADAM file
    val reads: RDD[ADAMRecord] = Input(sc, args.readInput, reference, config)

    // create stats/config item
    val stats = new AvocadoConfigAndStats(sc, args.debug, reads, reference)

    // apply read translation steps
    log.info("Processing reads.")
    val cleanedReads = preProcessReads(reads)

    // initial assignment of reads to variant calling algorithms
    log.info("Partitioning reads.")
    val partitionedReads = partitionReads(cleanedReads, stats)

    // call variants on filtered reads and pileups
    log.info("Calling variants.")
    val calledVariants = callVariants(partitionedReads)

    // post process variants
    log.info("Post-processing variants.")
    val processedGenotypes: RDD[ADAMGenotype] = postProcessVariants(calledVariants, stats).flatMap(variantContext => variantContext.genotypes)

    // save variants to output file
    log.info("Writing calls to disk.")
    processedGenotypes.adamSave(args.variantOutput,
      args.blockSize,
      args.pageSize,
      args.compressionCodec,
      args.disableDictionary)
  }
}
