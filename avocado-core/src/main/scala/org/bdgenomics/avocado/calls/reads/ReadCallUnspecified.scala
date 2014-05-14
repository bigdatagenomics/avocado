/*
 * Copyright (c) 2013. Regents of the University of California
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

package org.bdgenomics.avocado.calls.reads

import org.bdgenomics.adam.avro.ADAMRecord
import org.bdgenomics.adam.models.ADAMVariantContext
import org.bdgenomics.avocado.calls.VariantCallCompanion
import org.bdgenomics.avocado.partitioners.PartitionSet
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.{ SparkContext, Logging }
import org.apache.spark.rdd.RDD

object ReadCallUnspecified extends VariantCallCompanion {

  val callName = "ReadUnspecified"

  def apply(stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration,
            partitions: PartitionSet): ReadCallUnspecified = {

    new ReadCallUnspecified()
  }
}

/**
 * Abstract class for calling variants on reads.
 */
class ReadCallUnspecified extends ReadCall {

  val companion = ReadCallUnspecified

  /**
   * Empty calling method.
   */
  def call(pileupGroups: RDD[ADAMRecord]): RDD[ADAMVariantContext] = {
    throw new IllegalArgumentException(companion.callName + " is not callable.")
  }

  // Call is generic, so is not callable
  override def isCallable() = false
}
