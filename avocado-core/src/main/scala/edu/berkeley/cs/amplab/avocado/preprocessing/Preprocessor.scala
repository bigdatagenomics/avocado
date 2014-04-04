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

package org.bdgenomics.avocado.preprocessing

import org.bdgenomics.adam.avro.ADAMRecord
import org.apache.commons.configuration.HierarchicalConfiguration
import org.apache.spark.rdd.RDD

object Preprocessor {

  private val stages = List(MarkDuplicates,
    RecalibrateBaseQualities,
    SortReads,
    CoalesceReads,
    RealignIndels)

  def apply(rdd: RDD[ADAMRecord],
            stageName: String,
            stageAlgorithm: String,
            config: HierarchicalConfiguration): RDD[ADAMRecord] = {

    // get configuration for this stage
    val stageConfig = config.configurationAt(stageName)

    // find and run stage
    val stage = stages.find(_.stageName == stageAlgorithm)

    assert(stage.isDefined, "Could not find stage with name: " + stageName)
    stage.get.apply(rdd, stageConfig)
  }

}
