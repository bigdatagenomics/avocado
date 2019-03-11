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
import org.bdgenomics.formats.avro.AlignmentRecord
import scala.collection.mutable.ListBuffer

/**
 * Utilities to hard limit coverage at a single genomic site, across all
 * genomic sites.
 */
private[avocado] object HardLimiter extends Serializable {

  /**
   * Maps over an RDD to limit coverage at a single site.
   *
   * @note Assumes that the input RDD is sorted. This invariant is checked
   *   at runtime.
   *
   * @param rdd The RDD to limit coverage on.
   * @param maxCoverageAtSite The maximum per site coverage
   * @return Returns the RDD with coverage downsampled.
   */
  def apply[T](rdd: RDD[(AlignmentRecord, T)],
               maxCoverageAtSite: Int): RDD[(AlignmentRecord, T)] = {

    rdd.mapPartitions(limitPartition(_, maxCoverageAtSite))
  }

  /**
   * Maps over an iterator and limits coverage at a given site.
   *
   * @note Assumes that the input iterator is sorted. This invariant is checked
   *   at runtime.
   *
   * @param iter The iterator to limit coverage over.
   * @param maxCoverageAtSite The maximum per site coverage
   * @return Returns the iterator with coverage downsampled.
   */
  private[util] def limitPartition[T](iter: Iterator[(AlignmentRecord, T)],
                                      maxCoverageAtSite: Int): Iterator[(AlignmentRecord, T)] = {

    // allocate state variables
    var buffer = ListBuffer.empty[(AlignmentRecord, T)]
    var bufferSize = 0

    // loop and pop reads into the buffer
    iter.flatMap(read => {

      val (emitted, newBuffer, newBufferSize) = processRead(read,
        buffer,
        bufferSize,
        maxCoverageAtSite)

      // update buffer and buffer size
      buffer = newBuffer
      bufferSize = newBufferSize

      // "i am a great man, but i am not a good man"
      //
      // we need to flatten the buffer on the last read
      // there's no good way to do this in scala, hence this code
      // which is a bad way to do this in scala
      //
      // c'est la vie
      if (iter.hasNext) {
        emitted
      } else {
        emitted ++ buffer
      }
    })
  }

  /**
   * Adds a read in to the buffer, if space allows.
   *
   * @param read Read to add to buffer.
   * @param buffer The current buffer.
   * @param priorSize The size of the current buffer.
   * @param maxSize The maximum allowable buffer size.
   * @return Returns a tuple containing (reads removed from the buffer for
   *   emission, the read buffer, the buffer size).
   */
  private[util] def processRead[T](
    read: (AlignmentRecord, T),
    buffer: ListBuffer[(AlignmentRecord, T)],
    priorSize: Int,
    maxSize: Int): (ListBuffer[(AlignmentRecord, T)], ListBuffer[(AlignmentRecord, T)], Int) = {

    // validate the start/end pos of this new read
    val readStart = read._1.getStart
    buffer.lastOption.foreach(kv => {
      val (lastRead, _) = kv
      assert(lastRead.getStart <= readStart,
        "New read (%s) is before last read (%s).".format(read, lastRead))
      assert(lastRead.getReferenceName == read._1.getReferenceName)
    })

    // any read that ends before this new read starts can be flushed
    var bufferSize = priorSize
    val (flushed, kept) = buffer.partition(kv => {
      val (r, _) = kv

      val canFlush = r.getEnd <= readStart

      // we are popping one read from the buffer
      if (canFlush) {
        bufferSize -= 1
      }

      canFlush
    })

    // if our buffer is smaller than the max size, we can add the read
    if (bufferSize < maxSize) {
      bufferSize += 1
      kept += read
    }

    (flushed, kept, bufferSize)
  }
}
