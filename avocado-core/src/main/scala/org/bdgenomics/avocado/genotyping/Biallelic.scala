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
package org.bdgenomics.avocado.genotyping

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.avocado.models.Observation
import org.bdgenomics.formats.avro.AlignmentRecord

/**
 * Produces genotype likelihoods using a biallelic model.
 */
object Biallelic extends Serializable {

  /**
   * Transforms an RDD of reads into an RDD of per allele/per sample likelihoods.
   *
   * @param rdd RDD of reads to generate likelihoods from.
   * @return Returns an RDD of Observations, keyed by the (site, allele, and
   *   sample ID) that was observed.
   */
  def observe(rdd: RDD[AlignmentRecord]): RDD[((ReferenceRegion, String, String), Observation)] = {
    rdd.flatMap(observeRead)
      .reduceByKey(_.merge(_))
  }

  /**
   * From a single read, emits likelihood observations.
   *
   * Emits a likelihood for each allele seen in this read.
   *
   * @param read Read to observe.
   * @return Returns an iterable collection of observations, keyed by (site,
   *   allele, sample ID).
   */
  private[genotyping] def observeRead(read: AlignmentRecord): Iterable[((ReferenceRegion, String, String), Observation)] = {
    ???
  }
}
