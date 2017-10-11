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
import org.bdgenomics.avocado.genotyping.{ DiscoverVariants => DiscoverVariantsFromRDD }
import org.bdgenomics.utils.cli._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

object MergeDiscovered extends BDGCommandCompanion {
  val commandName = "mergeDiscovered"
  val commandDescription = "Merge variants discovered from reads of multiple samples"

  def apply(cmdLine: Array[String]) = {
    new MergeDiscovered(Args4j[MergeDiscoveredArgs](cmdLine))
  }
}

class MergeDiscoveredArgs extends Args4jBase {
  @Argument(required = true,
    metaVar = "INPUT",
    usage = "A globbed path of all discovered variants",
    index = 0)
  var inputPath: String = null
  @Argument(required = true,
    metaVar = "OUTPUT",
    usage = "Location to write the variants",
    index = 1)
  var outputPath: String = null
}

class MergeDiscovered(
    protected val args: MergeDiscoveredArgs) extends BDGSparkCommand[MergeDiscoveredArgs] {

  val companion = MergeDiscovered

  def run(sc: SparkContext) {

    // load variants, drop duplicates, save
    sc.loadVariants(args.inputPath)
      .transformDataset(_.dropDuplicates("start",
        "end",
        "contigName",
        "referenceAllele",
        "alternateAllele"))
      .saveAsParquet(args.outputPath)
  }
}
