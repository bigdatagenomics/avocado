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
import spark.{RDD,SparkContext}
import spark.SparkContext._
import org.apache.hadoop.mapreduce.Job
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, AvroReadSupport}
import org.kohsuke.args4j.{Option => option, Argument}
import org.apache.hadoop.mapreduce.Job
import edu.berkeley.cs.amplab.adam.predicates.LocusPredicate
import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup, ADAMRecord, ADAMVariant}
import edu.berkeley.cs.amplab.adam.commands.{AdamSparkCommand,AdamCommandCompanion,ParquetArgs,SparkArgs,ReadProcessor}
import edu.berkeley.cs.amplab.adam.util.{Args4j,Args4jBase}
import edu.berkeley.cs.amplab.avocado.calls.pileup.{PileupCall,PileupCallSimpleSNP}
import edu.berkeley.cs.amplab.avocado.filters.pileup.{PileupFilter,PileupFilterOnMismatch}
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._ 
import edu.berkeley.cs.amplab.avocado.calls.reads.{ReadCall}

object Avocado extends AdamCommandCompanion {

  val commandName = "Avocado"
  val commandDescription = "Call variants using the Avocado system."

  def apply (args: Array[String]) = {
    new Avocado (Args4j[AvocadoArgs](args))
  }
}

class AvocadoArgs extends Args4jBase with ParquetArgs with SparkArgs {
  @Argument (metaVar = "READS", required = true, usage = "ADAM read-oriented data", index = 0)
  var readInput: String = _

  @Argument (metaVar = "VARIANTS", required = true, usage = "ADAM variant output", index = 1)
  var variantOutput: String = _
}

class Avocado (protected val args: AvocadoArgs) extends AdamSparkCommand [AvocadoArgs] {
  
  // companion object to this class - needed for AdamCommand framework
  val companion = Avocado
    
  /* TODO: Add these steps in. Currently do not have any read based calls implemented,
   * and are waiting on ADAM transformation pipeline to be completed.
  def filterReads ()
  {
  }

  def processReads ()
  {
  }
  */

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
		    pileups: Map [PileupCall, RDD [ADAMPileup]]): RDD [ADAMVariant] = {
    
    if (!reads.isEmpty && !pileups.isEmpty) {
      val readCalledVariants = reads.map (kv => kv._1.call (kv._2)).reduce (_ ++ _)
      
      val pileupCalledVariants = pileups.map (kv => kv._1.call (kv._2)).reduce (_ ++ _)
      
      readCalledVariants ++ pileupCalledVariants
    } else if (!reads.isEmpty) {
      reads.map (kv => kv._1.call (kv._2)).reduce (_ ++ _)
    } else if (!pileups.isEmpty) {
      pileups.map (kv => kv._1.call (kv._2)).reduce (_ ++ _)
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
    
    // load in reads from ADAM file
    val reads: RDD[ADAMRecord] = sc.adamLoad (args.readInput, Some (classOf[LocusPredicate]))
    
    // TODO: add in read filtering stage
    //val readsToCall = reads.filter (v => false)

    // translate non-filtered reads into pileups
    val readProcessor = new ReadProcessor
    val pileups: RDD[ADAMPileup] = reads.flatMap ({readProcessor.readToPileups})

    // apply filters to pileups - these 
    val pileupsToCall = filterPileups (pileups)

    // call variants on filtered reads and pileups
    // TODO: currently, we don't have read filtering or calling implemented, so we pass an empty list to the read input
    val calledVariants = callVariants (Map[ReadCall, RDD [ADAMRecord]](), pileupsToCall)

    // TODO: clean up variant call filters and add filtering hook here

    // save variants to output file
    calledVariants.adamSave (args.variantOutput)
  }

}
