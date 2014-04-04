/*
 * Copyright (c) 2013-2014. Regents of the University of California
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

package org.bdgenomics.avocado.calls

import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.rdd.RDD
import org.apache.spark.{ SparkContext, Logging }
import org.bdgenomics.adam.avro.{ ADAMRecord, ADAMPileup, ADAMVariant, ADAMGenotype }
import org.bdgenomics.adam.converters.GenotypesToVariantsConverter
import org.bdgenomics.adam.models.ADAMVariantContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

trait VariantCallCompanion {

  val callName: String

  protected def apply(stats: AvocadoConfigAndStats, config: SubnodeConfiguration): VariantCall

  final def apply(stats: AvocadoConfigAndStats,
                  globalConfig: HierarchicalConfiguration,
                  callSetName: String): VariantCall = {
    val config: SubnodeConfiguration = globalConfig.configurationAt(callSetName)

    apply(stats, config)
  }

}

/**
 * Abstract class for calling variants on reads.
 */
abstract class VariantCall extends Serializable with Logging {

  val companion: VariantCallCompanion

  def process(rdd: RDD[ADAMRecord]): RDD[ADAMRecord] = {
    var readsToProcess = rdd.cache

    readsToProcess
  }

  def processAndCall(rdd: RDD[ADAMRecord]): RDD[ADAMVariantContext] = {
    val processedReads = process(rdd)
    call(rdd)
  }

  final def genotypesToVariantContext(genotypes: List[ADAMGenotype],
                                      samples: Int = 1): List[ADAMVariantContext] = {

    val grouped = genotypes.groupBy(
      g => (g.getVariant.getContig.getContigId, g.getVariant.getVariantAllele))
      .map(kv => {
        val (k, g) = kv

        ADAMVariantContext.buildFromGenotypes(g)
      })

    grouped.toList
  }

  def call(rdd: RDD[ADAMRecord]): RDD[ADAMVariantContext]

  def isReadCall(): Boolean
  def isPileupCall(): Boolean
  def isCallable(): Boolean
}

