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
import org.bdgenomics.adam.models.{ ReferenceMapping, ReferenceRegion }
import org.bdgenomics.adam.rdd.RegionJoin
import org.bdgenomics.adam.rich.ReferenceMappingContext._
import org.bdgenomics.avocado.algorithms.debrujin.KmerGraph
import org.bdgenomics.avocado.models.{ Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.{ AlignmentRecord, NucleotideContigFragment }

object ReassemblyExplorer extends ExplorerCompanion {

  val explorerName: String = "ReassemblyExplorer"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Explorer = {
    val kmerLength = config.getInt("kmerLength", 20)
    new ReassemblyExplorer(kmerLength, stats.reference)
  }

  implicit object ContigReferenceMapping extends ReferenceMapping[(Long, NucleotideContigFragment)] with Serializable {
    override def getReferenceName(value: (Long, NucleotideContigFragment)): String = value._2.getContig.getContigName.toString
    override def getReferenceRegion(value: (Long, NucleotideContigFragment)): ReferenceRegion = ReferenceRegion(value._2).get
  }
}

import ReassemblyExplorer._

class ReassemblyExplorer(kmerLength: Int,
                         reference: RDD[NucleotideContigFragment]) extends Explorer with Logging {

  val companion: ExplorerCompanion = ReassemblyExplorer

  def discoverRegion(rr: (Iterable[AlignmentRecord], NucleotideContigFragment)): Iterable[Observation] = {
    val (reads, reference) = rr

    // extract the coordinates from the reference
    val coordinates = ReferenceRegion(reference).get

    // assemble us up a region
    val graphs = KmerGraph(kmerLength,
      Seq((coordinates, reference.getFragmentSequence)),
      reads.toSeq)

    // turn the reassembly graph into observations
    graphs.flatMap(_.toObservations)
  }

  def discover(reads: RDD[AlignmentRecord]): RDD[Observation] = {
    // zip reference contig fragments with uuids, cache
    val refIds = reference.zipWithUniqueId()
      .map(vk => (vk._2, vk._1))
      .cache()

    // filter mapped reads, join with reference contigs, then extract contig ids
    // ultimately, this should use the merge-sort join, not the broadcast join
    // will upgrade when ADAM-534 merges.
    val joinWithId = RegionJoin.partitionAndJoin(reference.context,
      refIds,
      reads.filter(_.getReadMapped))
      .map(kv => {
        (kv._1._1, kv._2)
      })

    // merge together reads that are from the same contig fragment, then
    // join against reference contigs, and throw away id
    // ideally, we should optimize the partitioning; this is a to-do for later
    val joinReadsAndContigs = joinWithId.groupByKey()
      .join(refIds)
      .map(kv => kv._2)

    // we are done with the original reference contig rdd, so we can unpersist
    refIds.unpersist()

    // reassemble our regions into observations
    joinReadsAndContigs.flatMap(discoverRegion)
  }
}
