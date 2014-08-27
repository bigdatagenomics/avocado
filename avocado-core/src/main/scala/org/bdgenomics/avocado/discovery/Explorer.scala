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

import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.rdd.RDD
import org.bdgenomics.avocado.models.AlleleObservation
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.AlignmentRecord

trait ExplorerCompanion {

  val explorerName: String

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Explorer

  final def apply(stats: AvocadoConfigAndStats,
                  globalConfig: HierarchicalConfiguration,
                  explorerSetName: String): Explorer = {
    val config: SubnodeConfiguration = globalConfig.configurationAt(explorerSetName)

    apply(stats, config)
  }
}

trait Explorer extends Serializable {

  val companion: ExplorerCompanion

  def discover(reads: RDD[AlignmentRecord]): RDD[AlleleObservation]
}
