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

package edu.berkeley.cs.amplab.avocado.stats

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.ADAMNucleotideContigFragment

private[stats] object GetReferenceContigLengths {
  
  /**
   * From a rdd of reference contigs, collects their lengths in an array.
   *
   * @param rdd RDD of reference contigs.
   * @return List of contig lengths.
   */
  def apply (rdd: RDD[ADAMNucleotideContigFragment]): Map[Int, Long] = {
    rdd.map(c => (c.getContigId.toInt, c.getContigLength.toLong)).collect.toMap
  }

}
