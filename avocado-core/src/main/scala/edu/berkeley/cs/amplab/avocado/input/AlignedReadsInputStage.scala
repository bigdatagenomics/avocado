/*
 * Copyright (c) 2014. Regents of the University of California
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

package org.bdgenomics.avocado.input

import org.bdgenomics.adam.avro.{ ADAMRecord, ADAMNucleotideContigFragment }
import org.bdgenomics.adam.predicates.LocusPredicate
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.ADAMRDDFunctions
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD

private[input] object AlignedReadsInputStage extends InputStage {

  val stageName = "AlignedReads"

  /**
   * Sets up and loads data using this input stage.
   *
   * @param inputPath Path to input files.
   * @param config Configuration for this input stage.
   * @param reference RDD containing reference information.
   * @return Returns an RDD of ADAM reads.
   */
  def apply(sc: SparkContext,
            inputPath: String,
            config: SubnodeConfiguration,
            reference: RDD[ADAMNucleotideContigFragment]): RDD[ADAMRecord] = {

    println("Loading reads in from " + inputPath)

    sc.adamLoad(inputPath, Some(classOf[LocusPredicate]))
  }

}
