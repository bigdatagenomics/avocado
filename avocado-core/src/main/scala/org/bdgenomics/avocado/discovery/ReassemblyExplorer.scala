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
import org.apache.spark.SparkContext._
import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{ ReferenceMapping, ReferenceRegion, SequenceDictionary }
import org.bdgenomics.adam.rdd.ShuffleRegionJoin
import org.bdgenomics.adam.rich.ReferenceMappingContext._
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.algorithms.debrujin.KmerGraph
import org.bdgenomics.avocado.models.{ Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.{ AlignmentRecord, NucleotideContigFragment }

object ReassemblyExplorer extends ExplorerCompanion {

  val explorerName: String = "ReassemblyExplorer"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Explorer = {
    val kmerLength = config.getInt("kmerLength", 20)
    new ReassemblyExplorer(kmerLength, stats.reference, stats.sequenceDict, stats.contigLengths)
  }

  implicit object ContigReferenceMapping extends ReferenceMapping[(Long, NucleotideContigFragment)] with Serializable {
    override def getReferenceName(value: (Long, NucleotideContigFragment)): String = value._2.getContig.getContigName.toString
    override def getReferenceRegion(value: (Long, NucleotideContigFragment)): ReferenceRegion = ReferenceRegion(value._2).get
  }
}

import ReassemblyExplorer._

class ReassemblyExplorer(kmerLength: Int,
                         reference: RDD[NucleotideContigFragment],
                         sd: SequenceDictionary,
                         contigLengths: Map[String, Long]) extends Explorer with Logging {

  val totalAssembledReferenceLength = contigLengths.values.sum

  val companion: ExplorerCompanion = ReassemblyExplorer

  def discoverRegion(rr: (Iterable[AlignmentRecord], NucleotideContigFragment)): Iterable[Observation] = RegionDiscovery.time {
    val (reads, reference) = rr

    // extract the coordinates from the reference
    val coordinates = ReferenceRegion(reference).get

    // assemble us up a region
    val graphs = BuildingGraph.time {
      KmerGraph(kmerLength,
        Seq((coordinates, reference.getFragmentSequence.toUpperCase)),
        reads.toSeq)
    }

    // turn the reassembly graph into observations
    ObservingGraph.time {
      graphs.flatMap(_.toObservations)
    }
  }

  def discover(reads: RDD[AlignmentRecord]): RDD[Observation] = {
    val joinReadsAndContigs = JoiningReads.time {
      // zip reference contig fragments with uuids, cache
      val refIds = reference.zipWithUniqueId()
        .map(vk => (vk._2, vk._1))
        .cache()

      // filter mapped reads, join with reference contigs, then extract contig ids
      val joinWithId = ShuffleRegionJoin.partitionAndJoin(reference.context,
        refIds,
        reads.filter(_.getReadMapped),
        sd,
        totalAssembledReferenceLength / reads.partitions.size)
        .flatMap(kv => {
          val ((fragmentId, fragment), read) = kv

          val fStart = fragment.getFragmentStartPosition
          val fEnd = fStart + fragment.getFragmentSequence.length

          if (fStart <= read.getStart && fEnd >= read.getEnd) {
            Some((kv._1._1, kv._2))
          } else {
            None
          }
        })

      // merge together reads that are from the same contig fragment, then
      // join against reference contigs, and throw away id
      // ideally, we should optimize the partitioning; this is a to-do for later
      val jrdd = joinWithId.groupByKey()
        .join(refIds)
        .map(kv => kv._2)

      // we are done with the original reference contig rdd, so we can unpersist
      refIds.unpersist()

      jrdd
    }

    // reassemble our regions into observations
    joinReadsAndContigs.flatMap(discoverRegion)
  }
}
