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
package org.bdgenomics.avocado.input

import org.bdgenomics.formats.avro.{ AlignmentRecord, NucleotideContigFragment }
import org.apache.commons.configuration.HierarchicalConfiguration
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD

object Input {

  // all our input stages
  val stages = List(AlignedReadsInputStage, SnapInputStage)

  /**
   * Builds the input stage that corresponds to the given stage name, and returns the read data
   * that the stage provides. The input stage name to use is collected from the provided
   * configuration.
   *
   * @param sc A SparkContext
   * @param inputPath Path to input read data.
   * @param reference
   * @param config Configuration file containing the necessary data.
   * @return Returns an RDD of read data.
   */
  def apply(sc: SparkContext,
            inputPath: String,
            reference: RDD[NucleotideContigFragment],
            config: HierarchicalConfiguration): RDD[AlignmentRecord] = {
    // get input stage to use; if none is specified, default to input being aligned reads
    val stageName: String = config.getString("inputStage", "AlignedReads")

    val stage = stages.find(_.stageName == stageName)

    stage match {
      case Some(s: InputStage) => {
        val stageConfig = config.configurationAt(stageName)

        s.apply(sc, inputPath, stageConfig, reference)
      }
      case None => {
        throw new IllegalArgumentException("No input stage with name: " + stageName)
      }
    }
  }
}
