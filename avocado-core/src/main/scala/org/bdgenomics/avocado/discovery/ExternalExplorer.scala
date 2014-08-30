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
import org.bdgenomics.avocado.models.{ Observation, ReadObservation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.AlignmentRecord

object ExternalExplorer extends ExplorerCompanion {

  val explorerName: String = "ExternalExplorer"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Explorer = {
    new ExternalExplorer()
  }
}

class ExternalExplorer extends Explorer {

  val companion: ExplorerCompanion = ExternalExplorer

  def discover(reads: RDD[AlignmentRecord]): RDD[Observation] = {
    reads.map(r => new ReadObservation(r).asInstanceOf[Observation])
  }
}
