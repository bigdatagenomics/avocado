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

package edu.berkeley.cs.amplab.avocado.calls

import org.apache.commons.configuration.{HierarchicalConfiguration, SubnodeConfiguration}
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.{ADAMRecord, ADAMPileup, ADAMVariant, ADAMGenotype}
import edu.berkeley.cs.amplab.adam.models.ADAMVariantContext
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._
import edu.berkeley.cs.amplab.avocado.calls.pileup.{PileupCallSimpleSNP, PileupCallUnspecified}
import edu.berkeley.cs.amplab.avocado.calls.reads.{ReadCallAssemblyPhaser, ReadCallUnspecified}
import edu.berkeley.cs.amplab.avocado.stats.AvocadoConfigAndStats

object VariantCaller {

  private val calls = List(PileupCallSimpleSNP,
                           PileupCallUnspecified,
                           ReadCallAssemblyPhaser,
                           ReadCallUnspecified)

  def apply (callName: String,
             callAlgorithm: String, 
             stats: AvocadoConfigAndStats, 
             config: HierarchicalConfiguration): VariantCall = {
    val call = calls.find (_.callName == callAlgorithm)

    call match {
      case Some(c) => {
        c.apply(stats, config, callName)
      }
      case None => throw new IllegalArgumentException("Invalid variant calling algorithm provided.")
    }
  }

}
