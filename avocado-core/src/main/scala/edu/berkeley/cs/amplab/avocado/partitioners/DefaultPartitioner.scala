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

package edu.berkeley.cs.amplab.avocado.partitioners

import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.adam.models.ReferenceRegion
import edu.berkeley.cs.amplab.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import scala.collection.immutable.{SortedMap, TreeSet, SortedSet}

object DefaultPartitioner extends ReferencePartitionerCompanion {

  val partitionerName = "defaultPartitioner"

  def apply (rdd: RDD[ADAMRecord],
             subnodeConfiguration: SubnodeConfiguration,
             stats: AvocadoConfigAndStats): ReferencePartitioner = {

    new DefaultPartitioner(stats.contigLengths)
  }
}

class DefaultPartitioner (contigLengths: Map[Int, Long]) extends ReferencePartitioner {

  val companion = DefaultPartitioner

  def computePartitions (): PartitionSet = {
    var counter = 0
    var partitions = SortedMap[ReferenceRegion, Int]()
    
    contigLengths.map(kv => {
      val (id, len) = kv
      
      for (i <- 0 until (len / 1000).toInt) {
        partitions += (ReferenceRegion(id, i * 1000, i * 1000 + 999) -> counter)
        counter += 1
      }
    })

    new PartitionSet(partitions)
  }

}
