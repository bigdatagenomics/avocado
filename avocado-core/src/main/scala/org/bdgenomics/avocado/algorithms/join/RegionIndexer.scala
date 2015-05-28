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
package org.bdgenomics.avocado.algorithms.join

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion._
import org.bdgenomics.adam.models.{ ReferencePosition, ReferenceRegion }
import scala.annotation.tailrec

private[join] object RegionIndexer extends Serializable {

  def buildBySamplingRDD[T](rdd: RDD[(ReferenceRegion, T)],
                            fraction: Double = 0.01,
                            granularity: Long = 100000L): RegionIndexer = {
    val partitions = rdd.partitions.length

    val sampledCounts: Seq[(ReferencePosition, Long)] = rdd.sample(false, fraction)
      .map(kv => (ReferencePosition(kv._1.referenceName, (kv._1.start / granularity) * granularity), 1L))
      .reduceByKey(_ + _)
      .collect()
      .toSeq
      .sortBy(_._1)

    val countsPerBucket = (sampledCounts.map(_._2)
      .sum
      .toDouble / partitions.toDouble).toLong

    @tailrec def build(iter: Iterator[(ReferencePosition, Long)],
                       l: List[ReferenceRegion],
                       currCtg: String,
                       currStart: Long,
                       currEnd: Long,
                       currCount: Long): Seq[ReferenceRegion] = {
      if (!iter.hasNext) {
        (ReferenceRegion(currCtg, currStart, currEnd) :: l).toSeq
      } else {
        val (pos, count) = iter.next
        val end = pos.pos + granularity
        val updatedCount = currCount + count

        val (nextL, nextCtg, nextStart, nextEnd, nextCount) = if (pos.referenceName == currCtg) {
          (ReferenceRegion(currCtg, currStart, currEnd) :: l,
            pos.referenceName,
            pos.pos,
            pos.pos + granularity,
            updatedCount)
        } else {
          if (updatedCount > countsPerBucket) {
            val buckets = (updatedCount.toDouble / countsPerBucket.toDouble)
            val gap = (currEnd - currStart) / buckets

            ((0 until buckets.toInt).map(i => {
              ReferenceRegion(currCtg,
                currStart + (i * gap).toLong,
                currStart + ((i + 1) * gap).toLong)
            }).toList ::: l,
              currCtg,
              (currEnd - (buckets - buckets.toLong.toDouble) * gap).toLong,
              currEnd,
              (currCount - buckets * countsPerBucket).toLong)
          } else {
            (l,
              currCtg,
              currStart,
              pos.pos + granularity,
              updatedCount)
          }
        }

        build(iter, nextL, nextCtg, nextStart, nextEnd, nextCount)
      }
    }

    val iter = sampledCounts.toIterator
    val (firstPos, firstCount) = iter.next

    RegionIndexer(build(iter,
      List.empty,
      firstPos.referenceName,
      firstPos.pos,
      firstPos.pos + granularity,
      firstCount))
  }
}

private[join] case class RegionIndexer(regions: Seq[ReferenceRegion]) extends Indexer {

  val regionMap = regions.zipWithIndex
    .groupBy(kv => kv._1.referenceName)

  def numPartitions: Int = regions.length

  def getPartitions(region: ReferenceRegion): Iterable[Int] = {
    require(regionMap.contains(region.referenceName),
      "Couldn't find reference region for %s in (%s).".format(region,
        regionMap.keys.mkString(",")))

    regionMap(region.referenceName)
      .filter(kv => kv._1.overlaps(region))
      .map(_._2)
  }
}
