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

private[postprocessing] object FilterStrandBias extends PostprocessingStage {

  val stageName = "filterStrandBias"

  /**
   * Filters out genotype calls that show signs of strand bias.
   *
   * @param rdd Rdd on which to filter.
   * @param stats Global stats.
   * @param config Config from which to pull config data.
   * @return A filtered RDD.
   */
  def apply(rdd: RDD[VariantContext],
            stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): RDD[VariantContext] = {

    val pvalue = config.getFloat("fisherStrandBiasPValue")

    val genotypeFilter = new StrandBiasFilter(pvalue)

    genotypeFilter.filter(rdd)
  }
}

private[postprocessing] class StrandBiasFilter(strandBiasPValueThreshold: Float) extends GenotypeAttributeFilter[Float] {

  val filterName = "FISHER_STRAND_BIAS_PVALUE>=%f".format(strandBiasPValueThreshold)

  def keyFn(g: Genotype): Option[Float] = Option(g.getVariantCallingAnnotations)
    .flatMap(a => {
      val f: java.lang.Float = a.getFisherStrandBiasPValue
      Option(f).map(_.toFloat)
    })

  def filterFn(strandBiasPValue: Float): Boolean = strandBiasPValue >= strandBiasPValueThreshold
}
