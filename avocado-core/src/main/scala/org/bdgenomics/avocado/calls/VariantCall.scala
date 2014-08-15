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
package org.bdgenomics.avocado.calls

import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.rdd.RDD
import org.apache.spark.{ SparkContext, Logging }
import org.bdgenomics.formats.avro.{ AlignmentRecord, Pileup, Variant, Genotype }
import org.bdgenomics.adam.converters.GenotypesToVariantsConverter
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.partitioners.PartitionSet
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

trait VariantCallCompanion {

  val callName: String

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration,
                      partitions: PartitionSet): VariantCall

  final def apply(stats: AvocadoConfigAndStats,
                  globalConfig: HierarchicalConfiguration,
                  callSetName: String,
                  partitions: PartitionSet): VariantCall = {
    val config: SubnodeConfiguration = globalConfig.configurationAt(callSetName)

    apply(stats, config, partitions)
  }

}

/**
 * Abstract class for calling variants on reads.
 */
abstract class VariantCall extends Serializable with Logging {

  val companion: VariantCallCompanion

  def process(rdd: RDD[AlignmentRecord]): RDD[AlignmentRecord] = {
    var readsToProcess = rdd.cache

    readsToProcess
  }

  def processAndCall(rdd: RDD[AlignmentRecord]): RDD[VariantContext] = {
    val processedReads = process(rdd)
    call(rdd)
  }

  final def genotypesToVariantContext(genotypes: List[Genotype],
                                      samples: Int = 1): List[VariantContext] = {

    val grouped = genotypes.groupBy(
      g => (g.getVariant.getContig.getContigName.toString, g.getVariant.getAlternateAllele))
      .map(kv => {
        val (k, g) = kv

        VariantContext.buildFromGenotypes(g)
      })

    grouped.toList
  }

  def call(rdd: RDD[AlignmentRecord]): RDD[VariantContext]

  def isCallable(): Boolean
}

