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
import org.bdgenomics.adam.rdd.ADAMSaveAnyArgs
import org.bdgenomics.avocado.realigner.{ ConsensusRealigner, Realigner }
import org.bdgenomics.utils.cli._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

object Reassemble extends BDGCommandCompanion {
  val commandName = "reassemble"
  val commandDescription = "Reassemble reads to canonicalize variants"

  def apply(cmdLine: Array[String]) = {
    new Reassemble(Args4j[ReassembleArgs](cmdLine))
  }
}

class ReassembleArgs extends Args4jBase with ADAMSaveAnyArgs with ParquetArgs {
  @Argument(required = true,
    metaVar = "INPUT",
    usage = "The ADAM, BAM or SAM file to reassemble",
    index = 0)
  var inputPath: String = null
  @Argument(required = true,
    metaVar = "OUTPUT",
    usage = "Location to write the transformed data",
    index = 1)
  var outputPath: String = null
  @Argument(required = true,
    metaVar = "KMER_LENGTH",
    usage = "The k-mer length to use for reassembly",
    index = 2)
  var kmerLength = -1
  @Args4jOption(required = false,
    name = "-single",
    usage = "Saves OUTPUT as single file")
  var asSingleFile: Boolean = false
  @Args4jOption(required = false,
    name = "-defer_merging",
    usage = "Defers merging single file output")
  var deferMerging: Boolean = false
  @Args4jOption(required = false,
    name = "-use_consensus_realigner",
    usage = "If provided, uses consensus-mode realigner.")
  var useConsensusRealigner: Boolean = false

  // required by ADAMSaveAnyArgs
  var sortFastqOutput: Boolean = false
}

class Reassemble(
    protected val args: ReassembleArgs) extends BDGSparkCommand[ReassembleArgs] {

  val companion = Reassemble

  def run(sc: SparkContext) {

    // k-mer length must be positive
    require(args.kmerLength > 0,
      "k-mer length must be a positive value (%d given).".format(
        args.kmerLength))

    // load reads
    val reads = sc.loadAlignments(args.inputPath)

    // realign the reads
    val reassembledReads = if (args.useConsensusRealigner) {
      ConsensusRealigner.realign(reads, args.kmerLength)
    } else {
      Realigner.realign(reads, args.kmerLength)
    }

    // save the reads
    reassembledReads.save(args)
  }
}
