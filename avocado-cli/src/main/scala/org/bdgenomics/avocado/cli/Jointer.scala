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

import htsjdk.samtools.ValidationStringency
import org.apache.spark.SparkContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.ADAMSaveAnyArgs
import org.bdgenomics.avocado.genotyping.{
  JointAnnotatorCaller,
  SquareOffReferenceModel
}
import org.bdgenomics.utils.cli._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

object Jointer extends BDGCommandCompanion {
  val commandName = "jointer"
  val commandDescription = "Joint call and annotate variants."

  def apply(cmdLine: Array[String]) = {
    new Jointer(Args4j[JointerArgs](cmdLine))
  }
}

class JointerArgs extends Args4jBase with ADAMSaveAnyArgs with ParquetArgs {
  @Argument(required = true,
    metaVar = "INPUT",
    usage = "A globbed path of all genotyped variants",
    index = 0)
  var inputPath: String = null

  @Argument(required = true,
    metaVar = "OUTPUT",
    usage = "Location to write the squared off VCF",
    index = 1)
  var outputPath: String = null

  @Args4jOption(required = false,
    name = "-from_gvcf",
    usage = "Provide if input is gVCF with reference model, and not already squared off.")
  var fromGvcf: Boolean = false

  @Args4jOption(required = false,
    name = "-single",
    usage = "Save as a single VCF file.")
  var asSingleFile: Boolean = false

  @Args4jOption(required = false,
    name = "-defer_merging",
    usage = "Defers merging single file output.")
  var deferMerging: Boolean = false

  @Args4jOption(required = false,
    name = "-disable_fast_concat",
    usage = "Disables the parallel file concatenation engine.")
  var disableFastConcat: Boolean = false

  @Args4jOption(required = false,
    name = "-stringency",
    usage = "Stringency level for various checks; can be SILENT, LENIENT, or STRICT. Defaults to STRICT.")
  var stringency: String = "STRICT"

  // must be defined due to ADAMSaveAnyArgs, but unused here
  var sortFastqOutput: Boolean = false
}

class Jointer(
    protected val args: JointerArgs) extends BDGSparkCommand[JointerArgs] {

  val companion = Jointer

  def run(sc: SparkContext) {

    val stringency = ValidationStringency.valueOf(args.stringency)

    // load in input genotypes
    val genotypes = sc.loadGenotypes(args.inputPath)

    // are these in gVCF? if so, we must square off the variant matrix
    // once we've squared off, we can call
    if (args.fromGvcf) {

      val squaredOff = SquareOffReferenceModel(genotypes)

      // squaring off gives us sorted variants, so we need not resort at the end
      JointAnnotatorCaller(squaredOff)
        .saveAsVcf(args, stringency = stringency)
    } else {

      // load variants, drop duplicates, save
      JointAnnotatorCaller(genotypes)
        .transform(_.cache)
        .sort()
        .saveAsVcf(args, stringency = stringency)
    }
  }
}
