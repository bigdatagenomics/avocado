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

package edu.berkeley.cs.amplab.avocado.calls.reads

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.adam.models.ADAMVariantContext
import scala.collection.JavaConversions._
import org.scalatest.FunSuite
import scala.collection.mutable.{ArrayBuffer, Buffer, HashMap, HashSet, PriorityQueue, StringBuilder}

class ReadCallAssemblySuite extends FunSuite {
  
  test("Test the creation of several haplotype strings.") {
    val kms = ArrayBuffer[Kmer](new Kmer(new KmerPrefix("AC"), 'G', new KmerVertex, new KmerVertex),
                                new Kmer(new KmerPrefix("CG"), 'A', new KmerVertex, new KmerVertex),
                                new Kmer(new KmerPrefix("GA"), 'G', new KmerVertex, new KmerVertex))
    val kp = new KmerPath(kms)
    
    assert(kp.asHaplotypeString === "ACGAG")
  }

  test("Test log sum for similar values") {
    val hp = new HaplotypePair(null, null)
    val sum = hp.exactLogSumExp10(1.0, 1.0)
    
    assert(1.3 * 0.99 < sum && 1.3 * 1.01 > sum)
  }

  test("Test log sum for dissimilar values") {
    val hp = new HaplotypePair(null, null)
    val sum = hp.exactLogSumExp10(1.0, -3.0)
    
    assert(1.00004342 * 0.99 < sum && 1.00004342 * 1.01 > sum)
  }
}
