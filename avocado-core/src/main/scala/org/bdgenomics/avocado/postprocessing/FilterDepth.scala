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

private[postprocessing] object FilterDepth extends PostprocessingStage {

  val stageName = "filterDepth"

  /**
   * Filters out genotype calls that have low coverage. The coverage requirement can either
   * be specified as an absolute, or as a value relative to the measured coverage for this dataset.
   *
   * @throws IllegalArgumentException Throws exception if both config options are passed.
   *
   * @param rdd Rdd on which to filter.
   * @param stats Global stats.
   * @param config Config from which to pull config data.
   * @return A filtered RDD.
   */
  def apply(rdd: RDD[VariantContext],
            stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): RDD[VariantContext] = {

    val depth = if (config.containsKey("absoluteDepth") && !config.containsKey("relativeDepth")) {
      config.getInt("absoluteDepth")
    } else if (!config.containsKey("absoluteDepth") && config.containsKey("relativeDepth")) {
      (config.getDouble("relativeDepth") * stats.coverage).toInt
    } else if (config.containsKey("absoluteDepth") && config.containsKey("relativeDepth")) {
      throw new IllegalArgumentException("Both absoluteDepth and relativeDepth are specified.")
    } else {
      (0.75 * stats.coverage).toInt
    }

    val genotypeFilter = new DepthFilter(depth)

    genotypeFilter.filter(rdd)
  }
}

private[postprocessing] class DepthFilter(depth: Int) extends GenotypeFilter {

  /**
   * Filters genotypes that have low coverage.
   *
   * @param genotypes List of genotypes called at this site.
   * @return List of genotypes after filtering.
   */
  def filterGenotypes(genotypes: Seq[Genotype]): Seq[Genotype] = {
    val keyed = genotypes.map(g => (Option(g.getReadDepth), g))

    val genotypesNoStats: Seq[Genotype] = keyed.filter(t => t._1.isEmpty)
      .map(t => t._2)
    val genotypesWithStats: Seq[Genotype] = keyed.filter(t => t._1.isDefined)
      .filter(kv => kv._1.get >= depth)
      .map(kv => kv._2)

    genotypesNoStats ++ genotypesWithStats
  }

}
