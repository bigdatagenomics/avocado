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
package org.bdgenomics.avocado.discovery

import org.apache.commons.configuration.HierarchicalConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.models.Observation
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.AlignmentRecord

object Explore {

  val explorers: Seq[ExplorerCompanion] = Seq(ReadExplorer,
    ReassemblyExplorer,
    ExternalExplorer)

  def apply(explorerAlgorithm: String,
            explorerName: String,
            rdd: RDD[AlignmentRecord],
            stats: AvocadoConfigAndStats,
            config: HierarchicalConfiguration): RDD[Observation] = DiscoverObservations.time {

    explorers.find(_.explorerName == explorerAlgorithm)
      .fold(throw new IllegalArgumentException("Couldn't find exploration algorithm: " +
        explorerAlgorithm))(_.apply(stats,
        config,
        explorerName))
      .discover(rdd)
  }
}
