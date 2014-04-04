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

package org.bdgenomics.avocado.postprocessing

import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.avro.ADAMGenotype
import org.bdgenomics.adam.models.ADAMVariantContext
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

private[postprocessing] object FilterStrandBias extends PostprocessingStage {

  val stageName = "filterStrandBias"

  /**
   * Filters out genotype calls that have a strong strand bias. By default, filters
   * genotype calls that have strand bias greater than 0.75 or lower than 0.25. Parameters
   * can be set via config options "highLimit" and "lowLimit".
   *
   * @param rdd Rdd on which to filter.
   * @param stats Global stats.
   * @param config Config from which to pull config data.
   * @return A filtered RDD.
   */
  def apply(rdd: RDD[ADAMVariantContext],
            stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): RDD[ADAMVariantContext] = {

    val high = config.getDouble("highLimit", 0.75)
    val low = config.getDouble("lowLimit", 0.25)

    val genotypeFilter = new StrandBiasFilter(high, low)

    genotypeFilter.filter(rdd)
  }
}

private[postprocessing] class StrandBiasFilter(high: Double, low: Double) extends GenotypeFilter {

  /**
   * Filters genotypes that have a strand bias outside of acceptable range.
   *
   * @param genotypes List of genotypes called at this site.
   * @return List of genotypes after filtering.
   */
  def filterGenotypes(genotypes: Seq[ADAMGenotype]): Seq[ADAMGenotype] = {
    val keyed = genotypes.map(g => (Option(g.getReadDepth), Option(g.getReadsMappedForwardStrand), g))

    val genotypesNoStats: Seq[ADAMGenotype] = keyed.filter(t => t._1.isEmpty || t._2.isEmpty)
      .map(t => t._3)
    val genotypesWithStats: Seq[ADAMGenotype] = keyed.filter(t => t._1.isDefined && t._2.isDefined)
      .map(t => (t._2.get.toDouble / t._1.get.toDouble, t._3))
      .filter(kv => kv._1 > low && kv._1 < high)
      .map(kv => kv._2)

    genotypesNoStats ++ genotypesWithStats
  }
}
