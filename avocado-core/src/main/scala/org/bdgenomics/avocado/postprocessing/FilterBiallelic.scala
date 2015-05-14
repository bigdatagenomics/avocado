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
import org.bdgenomics.formats.avro.{ Genotype, VariantCallingAnnotations }
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

private[postprocessing] object FilterBiallelic extends PostprocessingStage {

  val stageName = "filterBiallelic"

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

    val genotypeFilter = new BiallelicFilter()

    genotypeFilter.filter(rdd)
  }
}

private[postprocessing] class BiallelicFilter extends GenotypeAttributeFilter[Boolean] {

  val filterName = "BIALLELIC"

  def keyFn(g: Genotype): Option[Boolean] = {
    Option(g.getSplitFromMultiAllelic)
  }

  def filterFn(isMultiAllelic: Boolean): Boolean = !isMultiAllelic
}
