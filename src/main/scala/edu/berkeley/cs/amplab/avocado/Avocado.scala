/*
 * Copyright (c) 2013. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.berkeley.cs.amplab.avocado

import parquet.hadoop.{ParquetOutputFormat, ParquetInputFormat}
import org.apache.spark.SparkContext._
import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.RDD
import org.apache.hadoop.mapreduce.Job
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, AvroReadSupport}
import org.kohsuke.args4j.{Option => option, Argument}
import org.apache.hadoop.mapreduce.Job
import edu.berkeley.cs.amplab.adam.predicates.LocusPredicate
import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup, ADAMRecord, ADAMVariant, ADAMGenotype}
import edu.berkeley.cs.amplab.adam.commands.{AdamSparkCommand, AdamCommandCompanion, ParquetArgs, SparkArgs}
import edu.berkeley.cs.amplab.adam.util.{Args4j, Args4jBase}
import edu.berkeley.cs.amplab.avocado.calls.pileup.{PileupCall, PileupCallSimpleSNP}
import edu.berkeley.cs.amplab.avocado.filters.pileup.{PileupFilter, PileupFilterOnMismatch}
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._ 
import edu.berkeley.cs.amplab.avocado.calls.reads.ReadCall

object Avocado extends AdamCommandCompanion {

  val commandName = "Avocado"
  val commandDescription = "Call variants using the Avocado system."

  var coverage = 1
  var reducePartitions = 1
  
  def apply (args: Array[String]) = {
    new Avocado (Args4j[AvocadoArgs](args))
  }
}

class AvocadoArgs extends Args4jBase with ParquetArgs with SparkArgs {
  @Argument (metaVar = "READS", required = true, usage = "ADAM read-oriented data", index = 0)
  var readInput: String = _

  @Argument (metaVar = "VARIANTS", required = true, usage = "ADAM variant output", index = 1)
  var variantOutputV: String = _

  @Argument (metaVar = "GENOTYPES", required = true, usage = "ADAM genotype output", index = 2)
  var variantOutputG: String = _
}

class Avocado (protected val args: AvocadoArgs) extends AdamSparkCommand [AvocadoArgs] with Logging {
  
  initLogging ()

  // companion object to this class - needed for AdamCommand framework
  val companion = Avocado
    
  /* TODO: Add these steps in. Currently do not have any read based calls implemented.
  def filterReads ()
  {
  }*/

  /**
   * Applies several pre-processing steps to the read pipeline. Currently, these are the default
   * steps in the ADAM processing pipeline.
   *
   * @param[in] reads RDD of reads to process.
   * @return RDD containing reads that have been sorted and deduped.
   */
  def processReads (reads: RDD[ADAMRecord]): RDD[ADAMRecord] = {
    
    // must sort reads before deduping them
    val sortedReads = reads.adamSortReadsByReferencePosition ()

    // dedupes reads across chromosomes
    val dedupedReads = sortedReads.adamMarkDuplicates ()

    // TODO: add in BQSR and indel realignment
    // BQSR: need VCF to mark out sites
    // Indel Realignment: waiting on algorithm in adam
    
    dedupedReads
  }

  /**
   * Function to pass over all pileups. Will either mark them for calling and provide a calling algorithm
   * for them, or will discard them.
   *
   * @param[in] pileups A Spark RDD containing pileups which are ready for calling.
   */
  def filterPileups (pileups: RDD[ADAMPileup]): Map [PileupCall, RDD [ADAMPileup]] = {
    
    // TODO: currently only supports a single filter and calling algorithm - need to extend
    val call = new PileupCallSimpleSNP
    val filter = new PileupFilterOnMismatch
    
    // run filter over input pileups, and assign a variant calling algorithm to them
    val map: Map [PileupCall, RDD [ADAMPileup]] = Map (call -> filter.filter (pileups))
    
    map.foreach(kv => log.info("For " + kv._1.callName + ", " + kv._2.count + " pileups to call."))
    
    map
  }

  /**
   * Applies variant calling algorithms to reads and pileups. Reduces down and returns called variants.
   *
   * @param[in] reads Map containing read calling algorithm and RDD record pairs.
   * @param[in] pileups Map containing pileup calling algorithm and RDD pileup pairs.
   * @return Joined output of variant calling algorithms.
   */
  def callVariants (reads: Map [ReadCall, RDD [ADAMRecord]],
		    pileups: Map [PileupCall, RDD [ADAMPileup]]): RDD [(ADAMVariant, List[ADAMGenotype])] = {
    
    if (!reads.isEmpty && !pileups.isEmpty) {
      log.info (reads.size.toString + " read calls.")
      log.info (pileups.size.toString + " pileup calls.")

      val readCalledVariants = reads.map (kv => {
        log.info ("Calling variants on reads using: " + kv._1.callName)

        kv._1.call (kv._2)
      }).reduce (_ ++ _)
      
      val pileupCalledVariants = pileups.map (kv => {
        log.info ("Calling variants on pileups using: " + kv._1.callName)

        kv._1.call (kv._2)
      }).reduce (_ ++ _)
      
      readCalledVariants ++ pileupCalledVariants
    } else if (!reads.isEmpty) {
      log.warn ("No pileups to call. Only calling reads.")
      log.info (reads.size.toString + " read calls.")

      reads.map (kv => {
        log.info ("Calling variants on reads using: " + kv._1.callName)

        kv._1.call (kv._2)
      }).reduce (_ ++ _)
    } else if (!pileups.isEmpty) {
      log.warn ("No reads to call. Only calling pileups.")
      log.info (pileups.size.toString + " pileup calls.")

      pileups.map (kv => {
        log.info ("Calling variants on pileups using: " + kv._1.callName)

        kv._1.call (kv._2)
      }).reduce (_ ++ _)
    } else {
      throw new IllegalArgumentException ("Must have some inputs to call.")
    }
  }

  /**
   * Main method. Implements body of variant caller. SparkContext and Hadoop Job are provided
   * by the AdamSparkCommand shell.
   *
   * @param[in] sc SparkContext for RDDs.
   * @param[in] job Hadoop Job container for file I/O.
   */
  def run (sc: SparkContext, job: Job) {

    log.info ("Starting avocado...")
    log.info ("Loading reads in from " + args.readInput)
    
    // load in reads from ADAM file
    val reads: RDD[ADAMRecord] = sc.adamLoad (args.readInput, Some (classOf[LocusPredicate]))
    val readCount = reads.count

    log.info ("Read " + readCount + " reads from " + args.readInput)

    // apply read translation steps
    log.info ("Processing reads.")
    val cleanedReads = processReads (reads)
    
    // TODO: add in read filtering stage
    //val readsToCall = reads.filter (v => false)

    log.info ("Translating reads into pileups.")

    // translate non-filtered reads into pileups
    val pileups: RDD[ADAMPileup] = cleanedReads.adamRecords2Pileup ()
    val pileupCount = pileups.count ()
    Avocado.coverage = (pileupCount / reads.count ()).toInt
    Avocado.reducePartitions = pileups.partitions.length * Avocado.coverage / 2

    log.info ("Have " + pileupCount + " pileups at " + Avocado.coverage + "x coverage. Aggregating...")

    val aggregatedPileups = pileups.adamAggregatePileups (Avocado.reducePartitions)
    val aggregatedCount = aggregatedPileups.count ()

    log.info ("Aggregation reduces pileups from " + pileupCount + " to " + aggregatedCount)

    log.info ("Applying filters to " + aggregatedCount + " pileups.")

    // apply filters to pileups - these 
    val pileupsToCall = filterPileups (aggregatedPileups)

    log.info ("Calling variants on pileups.")

    // call variants on filtered reads and pileups
    // TODO: currently, we don't have read filtering or calling implemented, so we pass an empty list to the read input
    val calledVariants = callVariants (Map[ReadCall, RDD [ADAMRecord]](), pileupsToCall)
    val variantCount = calledVariants.count ()

    // TODO: clean up variant call filters and add filtering hook here
    log.info ("Writing out " + variantCount + " variants.")

    // save variants to output file
    calledVariants.map (_._1).adamSave (args.variantOutputV)
    calledVariants.flatMap (_._2).adamSave (args.variantOutputG)
  }

}
