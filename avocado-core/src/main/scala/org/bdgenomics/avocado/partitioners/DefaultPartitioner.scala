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
package org.bdgenomics.avocado.partitioners

import org.bdgenomics.adam.avro.ADAMRecord
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import scala.collection.immutable.{ SortedMap, TreeSet, SortedSet }

/**
 * Configuration object for a default partitioner. The default partitioner performs no region
 * filtering, and just creates evenly sized buckets across all chromosomes. This partitioner
 * can be used for ensuring that all reads are mapped.
 *
 * This partitioner takes one optional parameter:
 * - bucketSize, Int: This parameter sets the size (in terms of reference positions) of a single
 * bucket of reads.
 */
object DefaultPartitioner extends ReferencePartitionerCompanion {

  val partitionerName = "defaultPartitioner"

  /**
   * Creates a default partitioner.
   *
   * @param rdd RDD to partition.
   * @param subnodeConfiguration Configuration for partitioner.
   * @param stats Statistics module for configuring partitioner.
   * @return Returns a default partitioner.
   */
  def apply(rdd: RDD[ADAMRecord],
            subnodeConfiguration: SubnodeConfiguration,
            stats: AvocadoConfigAndStats): ReferencePartitioner = {

    val bucketSize = subnodeConfiguration.getInt("bucketSize", 1000)

    new DefaultPartitioner(stats.contigLengths, bucketSize)
  }
}

/**
 * A default partitioner cuts the genome up into equally sized partitions which are spread
 * across reference contigs.
 *
 * @param contigLengths A map that maps a reference name to the length of that reference contig.
 * @param bucketSize The size (in reference positions) of a single bucket.
 */
class DefaultPartitioner(contigLengths: Map[String, Long],
                         bucketSize: Int) extends ReferencePartitioner {

  val companion = DefaultPartitioner

  /**
   * Computes the default partition set of equally sized buckets across all contigs.
   *
   * @return Partition set containing equally sized reference partitions.
   */
  def computePartitions(): PartitionSet = {
    var bucketCount: SortedMap[ReferenceRegion, Int] = SortedMap[ReferenceRegion, Int]()

    contigLengths.foreach(kv => {
      val (contigName, contigLength) = kv

      // if contig doesn't divide perfectly by bucket size, then add an extra bucket
      val buckets = if (contigLength % bucketSize == 0) {
        contigLength / bucketSize
      } else {
        contigLength / bucketSize + 1
      }

      bucketCount += (ReferenceRegion(contigName, 0L, contigLength) -> buckets.toInt)
    })

    // scan across map to get proper offsets - scan doesn't quite work correctly on maps...
    val offsets = bucketCount.values.scan(0)((a: Int, b: Int) => a + b).dropRight(1)
    val keys = bucketCount.keys
    val partitionOffsets = SortedMap[ReferenceRegion, Int]() ++ keys.zip(offsets)

    new DefaultPartitionSet(partitionOffsets, bucketSize)
  }

}

/**
 * Optimized PartitionSet implementation for default partitioner. Optimizes by computing
 * partition IDs arithmetically, which eliminates need to iterate across large map.
 *
 * @param mapping Mapping of contigs to the first bucket ID in that set.
 * @param bucketSize Size of a single bucket.
 */
class DefaultPartitionSet(mapping: SortedMap[ReferenceRegion, Int],
                          bucketSize: Int) extends PartitionSet(mapping) {

  // maps contig ID to the id of the first bucket in that contig
  val contigIdMap: Map[String, Int] = mapping.map(kv => (kv._1.referenceName, kv._2))

  // maps contig ID to the reference region of that contig
  val contigRegions: Map[String, ReferenceRegion] = mapping.map(kv => (kv._1.referenceName, kv._1))

  /**
   * Overrides main implementation with a faster implementation that does not check
   * the set of all partitions. Rather, it just computes the partition mapping through
   * arithmetic.
   *
   * @param region Reference region to find partition for.
   * @return IDs of all partitions the region maps to.
   */
  override def getPartition(region: ReferenceRegion): List[Int] = {
    try {
      // check that this region is inside of our reference contig
      if (contigRegions(region.referenceName).contains(region)) {
        val offset = contigIdMap(region.referenceName)
        var ids = List[Int]()

        // step across region by bucket size
        for (i <- region.start to region.end by bucketSize.toLong) {
          val bucket = i.toInt / bucketSize + offset

          ids = bucket :: ids
        }

        // if end of region is in another region that has not already been added, add
        val endBucket = region.end.toInt / bucketSize + offset
        if (!ids.contains(endBucket)) {
          ids = endBucket :: ids
        }

        ids
      } else {
        throw new IllegalArgumentException("Region on contig " + region.referenceName + " with start at " + region.start +
          " and end at " + region.end + " is outside of the reference contig.")
      }
    } catch {
      case _: Throwable => throw new IllegalArgumentException("Received region with contig ID " + region.referenceName +
        ", but do not have a matching reference contig.")
    }
  }

  /**
   * Reads are always inside of the reference.
   *
   * @param region Reference region to check.
   * @return Always returns true.
   */
  override def isInSet(region: ReferenceRegion): Boolean = true

  /**
   * Reads are never outside of the reference.
   *
   * @param region Reference region to check.
   * @return Always returns true.
   */
  override def isOutsideOfSet(region: ReferenceRegion): Boolean = false
}
