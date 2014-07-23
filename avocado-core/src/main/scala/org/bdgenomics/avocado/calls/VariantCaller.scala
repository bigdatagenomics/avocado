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
import org.bdgenomics.formats.avro.{ ADAMRecord, ADAMPileup, ADAMVariant, ADAMGenotype }
import org.bdgenomics.adam.models.ADAMVariantContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.calls.pileup.{ MPileupCallSimpleSNP, PileupCallSimpleSNP, PileupCallUnspecified }
import org.bdgenomics.avocado.calls.reads.{ ReadCallAssemblyPhaser, ReadCallUnspecified, ReadCallExternal }
import org.bdgenomics.avocado.partitioners.PartitionSet
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats

object VariantCaller {

  private val calls = List(PileupCallSimpleSNP,
    MPileupCallSimpleSNP,
    PileupCallUnspecified,
    ReadCallAssemblyPhaser,
    ReadCallUnspecified,
    ReadCallExternal)

  def apply(callName: String,
            callAlgorithm: String,
            stats: AvocadoConfigAndStats,
            config: HierarchicalConfiguration,
            partitions: PartitionSet): VariantCall = {
    val call = calls.find(_.callName == callAlgorithm)

    call match {
      case Some(c) => {
        c.apply(stats, config, callName, partitions)
      }
      case None => throw new IllegalArgumentException("Invalid variant calling algorithm provided.")
    }
  }

}
