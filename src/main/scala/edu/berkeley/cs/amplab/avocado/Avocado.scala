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
import edu.berkeley.cs.amplab.avocado.calls.pileup.{PileupCall, PileupCallSimpleSNP, PileupCallSNPVCFForMAF, PileupCallUnspecified}
import edu.berkeley.cs.amplab.avocado.filters.pileup.{PileupFilter, PileupFilterOnMismatch}
import edu.berkeley.cs.amplab.avocado.filters.reads.{ReadFilter, ReadFilterOnComplexity}
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._ 
import edu.berkeley.cs.amplab.avocado.calls.reads.{ReadCall, ReadCallUnspecified}
import edu.berkeley.cs.amplab.avocado.calls.VariantCall
import java.io.File

object Avocado extends AdamCommandCompanion {

  val commandName = "Avocado"
  val commandDescription = "Call variants using the Avocado system."

  var coverage = 15
  var reducePartitions = 15
  
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

  @option (name = "-m", usage = "VCF formatted file containing MAFs as GMAF attribute")
  var mafFileName: String = ""

  @option (name = "-bqsrVcf", usage = "VCF file for BQSR. If not provided, BQSR is not run, unless a VCF for MAF is provided.")
  var bqsrFile: File = null

  @option (name = "-locallyAssembleHighComplexity", usage = "Use local assembly for high complexity regions.")
  var locallyAssembleHighComplexity = false

  @option (name = "-locallyAssembleAll", usage = "Use local assembly for all reads. Takes precident over high complexity.")
  var locallyAssembleAll = false

  @option (name = "-aggregatePileups", usage = "Aggregates pileups wherever possible.")
  var aggregatePileups = false

  @option (name = "-coverageTolerance", usage = "Coverage tolerance for complexity filter. Default value = 0.5 (50%)")
  var coverageTolerance = 0.5
  
  @option (name = "-mappingThreshold", usage = "Mapping quality threshold for complexity filter. Default value = PHRED 30.")
  var mappingThreshold = 30
}

class Avocado (protected val args: AvocadoArgs) extends AdamSparkCommand [AvocadoArgs] with Logging {
  
  initLogging ()

  // companion object to this class - needed for AdamCommand framework
  val companion = Avocado
    
  /**
   * Assigns reads to a variant calling algorithm.
   *
   * @param reads An RDD of raw reads.
   * @return A map that maps RDDs of reads to variant calling algorithms.
   */
  def filterReads (reads: RDD[ADAMRecord]): Map[VariantCall, RDD[ADAMRecord]] = {
    if (args.locallyAssembleAll) {
      Map[VariantCall, RDD[ADAMRecord]](new ReadCallUnspecified -> reads)
    } else if (args.locallyAssembleHighComplexity) {
      // set coverage threshold
      val coverageHigh = (Avocado.coverage * (1.0 + args.coverageTolerance)).toInt
      val coverageLow = (Avocado.coverage * (1.0 - args.coverageTolerance)).toInt

      val filter = new ReadFilterOnComplexity (coverageHigh, coverageLow, args.mappingThreshold)
      filter.filter (reads)
    } else {
      // no read based calls
      Map[VariantCall, RDD[ADAMRecord]](new PileupCallUnspecified -> reads)
    }
  }

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
    
    // if vcf is provided for bqsr or maf, run bqsr - prefer bqsr argument
    val bqsrReads = if (args.bqsrFile != null) {
      dedupedReads.adamBQSR (args.bqsrFile)
    } else if (args.mafFileName != "") {
      dedupedReads.adamBQSR (new File (args.mafFileName))
    } else {
      dedupedReads
    }

    // TODO: add in indel realignment
    // Indel Realignment: waiting on algorithm in adam
    
    bqsrReads
  }

  /**
   * Function to pass over all pileups. Will either mark them for calling and provide a calling algorithm
   * for them, or will discard them.
   *
   * @param[in] pileups A Spark RDD containing pileups which are ready for calling.
   */
  def filterPileups (pileupMap: Map[PileupCall,RDD[ADAMPileup]]): Map [PileupCall, RDD [ADAMPileup]] = {
    // loop over callset
    pileupMap.map(kv => {
      val (call, pileups) = kv

      if (call.isCallable) {
        // if already assigned a viable pileup call, return
        kv
      } else {
        // if not, then provide a new call - if -m use PileupCallVCF
        val call = if (args.mafFileName == "") {
          new PileupCallSimpleSNP
        } else {
          new PileupCallSNPVCFForMAF(args.mafFileName)
        }
        
        val filter = new PileupFilterOnMismatch
        
        // run filter over input pileups, and assign a variant calling algorithm to them
        val filteredPileups: RDD [ADAMPileup] = filter.filter (pileups)
        
        (call, filteredPileups)
      }
    })
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
    
    // filter out generic calls that are not callable
    val readsFiltered = reads.filter (_._1.isCallable)
    val pileupsFiltered = pileups.filter (_._1.isCallable)

    if (!readsFiltered.isEmpty && !pileupsFiltered.isEmpty) {
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
    } else if (!readsFiltered.isEmpty) {
      log.warn ("No pileups to call. Only calling reads.")
      log.info (reads.size.toString + " read calls.")

      reads.map (kv => {
        log.info ("Calling variants on reads using: " + kv._1.callName)

        kv._1.call (kv._2)
      }).reduce (_ ++ _)
    } else if (!pileupsFiltered.isEmpty) {
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
   * Loops across all segregated calls. If call is a pileup based call, it converts its read set into pileups.
   *
   * @param reads Map that maps variant calling algorithms to RDDs of reads.
   * @return Tuple containing (Map of read based calls to RDDs of reads, Map of pileup based calls to RDDs of pileups).
   */
  def translateReadsToPileups (reads: Map[VariantCall,RDD[ADAMRecord]]): (Map[ReadCall,RDD[ADAMRecord]], Map[PileupCall,RDD[ADAMPileup]]) = {
    log.info ("Translating reads into pileups.")

    // separate out read and pileup calls
    val readCalls: Map[ReadCall,RDD[ADAMRecord]] = reads.filter(_._1.isReadCall)
      .map(kv => (kv._1.asInstanceOf[ReadCall], kv._2))
    val pileupCallsPreTranslation: Map[PileupCall,RDD[ADAMRecord]] = reads.filter(_._1.isPileupCall)
      .map(kv => (kv._1.asInstanceOf[PileupCall], kv._2))

    // translate non-filtered reads into pileups
    var pileupCalls: Map[PileupCall,RDD[ADAMPileup]] = pileupCallsPreTranslation.map(kv => {
      val (call, reads) = kv

      log.info ("Performing read to pileup conversion for " + call.callName + ".")

      val pileups: RDD[ADAMPileup] = reads.adamRecords2Pileup ()
      val pileupCount = pileups.count ()
            
      log.info ("Have " + pileupCount + " pileups at " + Avocado.coverage + "x coverage.")
      
      // aggregate pileups if desired
      val aggregatedPileups = if (args.aggregatePileups) {
        log.info ("Aggregating...")
        
        pileups.adamAggregatePileups (Avocado.reducePartitions)
      } else {
        pileups
      }

      (call, aggregatedPileups)
    })

    (readCalls, pileupCalls)
  }

  /**
   * Collects coverage statistics.
   *
   * @param reads RDD of reads.
   * @return Coverage derived from read data.
   */
  def getCoverage (reads: RDD[ADAMRecord]): Int = {
    // TODO: this is not accurate for exomes...
    val genomeStart: Long = reads.map(_.getStart.toLong).reduce(_ min _)
    val genomeEnd: Long = reads.map(_.end.get).reduce(_ max _)
    val baseCount = reads.map(_.getSequence.length.toLong).reduce(_ + _)

    val coverage = baseCount / (genomeEnd - genomeStart)

    coverage.toInt
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
    Avocado.coverage = getCoverage (reads)
    log.info ("Read " + readCount + " reads from " + args.readInput + " at " + Avocado.coverage + "x coverage.")

    // apply read translation steps
    log.info ("Processing reads.")
    val cleanedReads = processReads (reads)
    
    // initial assignment of reads to variant calling algorithms
    val filteredReads = filterReads (cleanedReads)

    // translate reads to pileups for pileup based algorithms
    val (readsToCall, pileupCalls) = translateReadsToPileups (filteredReads)

    // apply filters to pileups 
    val pileupsToCall = filterPileups (pileupCalls)

    // call variants on filtered reads and pileups
    val calledVariants = callVariants (readsToCall, pileupsToCall)
    val variantCount = calledVariants.count ()

    // TODO: clean up variant call filters and add filtering hook here
    log.info ("Writing out " + variantCount + " variants.")

    // save variants to output file
    calledVariants.map (_._1).adamSave (args.variantOutputV)
    calledVariants.flatMap (_._2).adamSave (args.variantOutputG)
  }

}
