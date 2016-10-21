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

import java.io.FileNotFoundException
import org.apache.hadoop.fs.Path
import org.apache.spark.SparkContext
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.formats.avro.AlignmentRecord
import scala.io.Source

private[avocado] object GapFilter {

  def apply(reads: AlignmentRecordRDD,
            gapPath: String): AlignmentRecordRDD = FilteringByGaps.time {

    // load gaps from file
    val gaps = loadGaps(reads.rdd.context,
      gapPath)

    // filter reads by gaps
    reads.transform(rdd => rdd.filter(filterGap(_, gaps)))
  }

  private[util] def loadGaps(sc: SparkContext,
                             filePath: String): Map[String, Seq[(Long, Long)]] = {

    // get underlying file system
    val path = new Path(filePath)
    val fs = Option(path.getFileSystem(sc.hadoopConfiguration))
      .getOrElse({
        throw new FileNotFoundException(
          "Couldn't find filesystem for %s with Hadoop configuration %s".format(
            path.toUri, sc.hadoopConfiguration))
      })

    // open stream from fs and make source
    val stream = fs.open(path)
    val source = Source.fromInputStream(stream)

    // parse lines
    val lines = source.getLines
      .map(l => {
        val cols = l.split("\\s+")

        require(cols.size == 3,
          "Expected 3 whitespace delimited columns in line %s.".format(l))

        (cols(0), cols(1).toLong, cols(2).toLong)
      }).toSeq

    // group lines by contig name, drop key, and return
    lines.groupBy(_._1)
      .map(kv => {
        val (contigName, tuples) = kv

        (contigName, tuples.map(t => (t._2, t._3)))
      })
  }

  private[util] def filterGap(read: AlignmentRecord,
                              gaps: Map[String, Seq[(Long, Long)]]): Boolean = {

    // get the gaps for the chromosome this read is aligned to
    val gapsOnMappedContig = gaps(read.getContigName)

    // does this read overlap any of the gaps?
    val gapsOverlapping = gapsOnMappedContig.filter(gap => {
      val (gapStart, gapEnd) = gap

      (gapStart < read.getEnd && gapEnd > read.getStart)
    })

    // discard this read if it overlapped a gap
    gapsOverlapping.isEmpty
  }
}

