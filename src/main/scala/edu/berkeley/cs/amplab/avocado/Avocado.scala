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
  
  val companion = Avocado
    
  /*def filterReads ()
  {
  }

  def processReads ()
  {
  }*/

  def filterPileups (pileups: RDD[ADAMPileup]): Map [PileupCall, RDD [ADAMPileup]] = {
    
    val call = new PileupCallSimpleSNP
    val filter = new PileupFilterOnMismatch

    val map: Map [PileupCall, RDD [ADAMPileup]] = Map (call -> filter.filter (pileups))
    
    map
  }

  def callVariants (reads: Map [ReadCall, RDD [ADAMRecord]],
		    pileups: Map [PileupCall, RDD [ADAMPileup]]): RDD [ADAMVariant] = {
  
    val readCalledVariants = reads.map (kv => kv._1.call (kv._2)).reduce (_ ++ _)

    val pileupCalledVariants = pileups.map (kv => kv._1.call (kv._2)).reduce (_ ++ _)

    readCalledVariants ++ pileupCalledVariants
  }

  def run (sc: SparkContext, job: Job) {
    
    val reads: RDD[ADAMRecord] = sc.adamLoad (args.readInput, Some (classOf[LocusPredicate]))
    
    // placeholder
    //val readsToCall = reads.filter (v => false)

    val readProcessor = new ReadProcessor
    val pileups: RDD[ADAMPileup] = reads.flatMap ({readProcessor.readToPileups})

    val pileupsToCall = filterPileups (pileups)

    // fixme: null
    val calledVariants = callVariants (null, pileupsToCall)

    println (calledVariants.count)
  }

}
