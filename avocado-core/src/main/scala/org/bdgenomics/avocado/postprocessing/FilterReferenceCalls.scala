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
import org.bdgenomics.formats.avro.{ Genotype, GenotypeAllele }
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

private[postprocessing] object FilterReferenceCalls extends PostprocessingStage {

  val stageName = "filterReferenceCalls"

  /**
   * Filters out genotype calls that are homozygous reference. If this filter is not used,
   * we emit a gVCF.
   *
   * @param rdd Rdd on which to filter.
   * @param stats Global stats.
   * @param config Config from which to pull config data.
   * @return A filtered RDD.
   */
  def apply(rdd: RDD[VariantContext],
            stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): RDD[VariantContext] = {

    val genotypeFilter = new ReferenceCallFilter()

    genotypeFilter.filter(rdd)
  }
}

private[postprocessing] class ReferenceCallFilter() extends GenotypeFilter {

  /**
   * Filters genotypes that are homozygous reference
   *
   * @param genotypes List of genotypes called at this site.
   * @return List of genotypes after filtering.
   */
  def filterGenotypes(genotypes: Seq[Genotype]): Seq[Genotype] = {
    val hasNonRefCall = genotypes.map(g => g.getAlleles
      .toArray
      .count(a => a == GenotypeAllele.Alt || a == GenotypeAllele.OtherAlt) > 0)
      .fold(false)(_ || _)

    if (hasNonRefCall) {
      genotypes
    } else {
      Seq()
    }
  }
}
