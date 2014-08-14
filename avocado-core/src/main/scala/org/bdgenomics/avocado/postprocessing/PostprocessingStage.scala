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
package org.bdgenomics.avocado.postprocessing

import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.Genotype
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

private[postprocessing] trait PostprocessingStage {

  val stageName: String

  def apply(rdd: RDD[VariantContext],
            stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): RDD[VariantContext]

}

private[postprocessing] trait GenotypeFilter extends Serializable {

  /**
   * Abstract method that must be implemented. Implements basic filtering on genotypes that
   * are inside a single variant context.
   *
   * @param genotypes Genotypes to filter.
   * @return Filtered genotypes.
   */
  def filterGenotypes(genotypes: Seq[Genotype]): Seq[Genotype]

  /**
   * Applies filtering and creates a new variant context, if called genotypes still exist.
   * If all genotypes have been filtered out, then an empty option (None) is returned.
   *
   * @param vc Variant context on which to filter.
   * @return If not all genotypes have been filtered out, a new variant context, else none.
   */
  def createNewVC(vc: VariantContext): Option[VariantContext] = {
    val filteredGt = filterGenotypes(vc.genotypes.toSeq)

    if (filteredGt.length > 0) {
      Some(VariantContext.buildFromGenotypes(filteredGt))
    } else {
      None
    }
  }

  /**
   * Applies the filtering described above across a full RDD.
   *
   * @param rdd RDD of variant contexts.
   * @return An RDD containing variant contexts after filtering.
   */
  def filter(rdd: RDD[VariantContext]): RDD[VariantContext] = {
    rdd.flatMap(vc => createNewVC(vc))
  }
}
