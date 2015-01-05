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
package org.bdgenomics.avocado.postprocessing

import org.apache.commons.configuration.HierarchicalConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

object Postprocessor {

  private val stages = List[PostprocessingStage](FilterDepth, FilterReferenceCalls)

  assert(stages.map(_.stageName).length == stages.map(_.stageName).distinct.length,
    "Postprocessing stages have duplicated names.")

  def apply(rdd: RDD[VariantContext],
            stageName: String,
            stageAlgorithm: String,
            stats: AvocadoConfigAndStats,
            config: HierarchicalConfiguration): RDD[VariantContext] = {

    val stage = stages.find(_.stageName == stageAlgorithm)

    stage match {
      case Some(s) => {
        val c = config.configurationAt(stageName)

        s.apply(rdd, stats, c)
      }
      case None => throw new IllegalArgumentException("Postprocessing stage " + stageAlgorithm + "does not exist.")
    }
  }

}
