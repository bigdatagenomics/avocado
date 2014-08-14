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
package org.bdgenomics.avocado.stats

import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichAlignmentRecord
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD

private[stats] object ScoreCoverage {

  // TODO: correct for whole genome without gaps, incorrect if gaps are present
  def apply(rdd: RDD[AlignmentRecord]): Double = {

    def coverageReducer(t1: (Long, Long, Long), t2: (Long, Long, Long)): (Long, Long, Long) = {
      (t1._1 min t2._1,
        t1._2 max t2._1,
        t1._3 + t2._3)
    }

    def readToParams(r: AlignmentRecord): (Long, Long, Long) = {
      val s = r.getStart
      val e = r.getEnd
      (s, e, e - s + 1)
    }

    val (start, end, bases) = rdd.filter(_.getReadMapped).map(readToParams)
      .reduce(coverageReducer(_, _))

    bases.toDouble / (end - start).toDouble
  }

}
