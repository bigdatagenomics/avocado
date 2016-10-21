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
package org.bdgenomics.avocado.util

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.GenomeBins
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.utils.misc.Logging
import scala.util.Random

private[avocado] object Downsampler extends Serializable with Logging {

  /**
   * Downsamples the number of reads in a given RDD partition.
   *
   * @param rdd RDD to downsample.
   * @param targetPartitionSize The maximum number of reads to have in a
   *   partition.
   * @param bins The genome bin descriptor.
   * @param optSeed An optional random seed.
   * @return Returns a downsampled RDD.
   */
  def downsample[T](
    rdd: RDD[(AlignmentRecord, T)],
    targetPartitionSize: Int,
    bins: GenomeBins,
    optSeed: Option[Int] = None): RDD[(AlignmentRecord, T)] = DownsampleReads.time {
    rdd.mapPartitionsWithIndex(downsamplePartition(_, _,
      targetPartitionSize.toDouble,
      bins,
      optSeed))
  }

  private def downsamplePartition[T](
    idx: Int,
    iter: Iterator[(AlignmentRecord, T)],
    targetPartitionSize: Double,
    bins: GenomeBins,
    optSeed: Option[Int]): Iterator[(AlignmentRecord, T)] = DownsamplePartition.time {

    val indicesToSkip = Set[Int]()
    /*4152, 2373, 290, 332, 4510, 4953,
       4511, 2346, 2537, 2670, 3247,
       3629, 289, 676)*/

    // get the region covered by this bin
    val region = bins.invert(idx)

    if (indicesToSkip(idx)) {

      log.warn("Skipping partition %d, which covers region %s.".format(idx, region))
      Iterator()
    } else {

      // get random number generator
      val rv = optSeed.fold(new Random)(v => new Random(v + idx))

      // how wide is this region?
      val binWidth = region.length.toDouble

      // how many reads have we seen?
      var readsSeen: Double = 0.0
      var readsEmitted: Double = 0.0

      def downsamplingRatio(read: AlignmentRecord): Option[Double] = {

        // how far are we into the region?
        val ratioThrough = (read.getStart - region.start).toDouble / binWidth

        // how many reads have we seen relative to what we expect to see?
        val ratioReadsSeen = readsSeen / targetPartitionSize

        Some(ratioReadsSeen / ratioThrough).filter(_ < 1.0)
      }

      val emittedGenotypes = iter.flatMap(readTuple => {

        // increment the read count
        readsSeen += 1.0

        // get the downsampling ratio
        val downsampleBy = downsamplingRatio(readTuple._1)

        // should we downsample or keep the read?
        val downsampleRead = downsampleBy.fold(false)(pctToKeep => {

          // sample the rv
          val sample = rv.nextDouble()

          sample > pctToKeep
        })

        if (downsampleRead) {
          None
        } else {
          readsEmitted += 1.0
          Some(readTuple)
        }
      })

      log.info("Saw %d reads for partition %d, emitted %d.".format(
        readsSeen.toInt, idx, readsEmitted.toInt))
      emittedGenotypes
    }
  }
}

