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

package edu.berkeley.cs.amplab.avocado.calls

import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup, ADAMVariant, ADAMGenotype}
import edu.berkeley.cs.amplab.adam.models.ADAMVariantContext
import edu.berkeley.cs.amplab.adam.projections.ADAMVariantField
import edu.berkeley.cs.amplab.adam.rdd.GenotypesToVariantsConverter
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._
import edu.berkeley.cs.amplab.avocado.Avocado
import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.RDD

/**
 * Abstract class for calling variants on reads. 
 */
abstract class VariantCall extends Serializable with Logging {

  val callName: String

  final def genotypesToVariantContext(genotypes: List[ADAMGenotype],
                                      samples: Int = 1): List[ADAMVariantContext] = {
    
    val conv = new GenotypesToVariantsConverter(false, false)

    val grouped = genotypes.groupBy(g => (g.getReferenceId, g.getAllele))
      .flatMap(kv => {
        val (k, g) = kv
        val sk = (k._1.toInt, k._2.toString)

        try {
          val v = conv.convertGenotypes(g, sk, None, Set[ADAMVariantField.Value](), g.length, samples)

          Some(new ADAMVariantContext(g.head.getPosition, List(v), g, None))
        } catch {
          case _ => None
        }
      })
    
    grouped.toList
  }

  def isReadCall (): Boolean
  def isPileupCall (): Boolean
  def isCallable (): Boolean
}

