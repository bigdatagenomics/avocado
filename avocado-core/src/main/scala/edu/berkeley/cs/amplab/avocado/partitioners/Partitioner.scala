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

package org.bdgenomics.avocado.partitioners

import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.avro.ADAMRecord
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

object Partitioner {

  val partitioners = List(DefaultPartitioner)

  def apply(rdd: RDD[ADAMRecord],
            globalConfig: HierarchicalConfiguration,
            partitionName: String,
            partitionerAlgorithm: String,
            stats: AvocadoConfigAndStats): ReferencePartitioner = {
    val partitioner = partitioners.find(_.partitionerName == partitionerAlgorithm)

    assert(partitioner.isDefined, "Could not find partitioner: " + partitionerAlgorithm)
    partitioner.get.apply(rdd, globalConfig, partitionName, stats)
  }

}
