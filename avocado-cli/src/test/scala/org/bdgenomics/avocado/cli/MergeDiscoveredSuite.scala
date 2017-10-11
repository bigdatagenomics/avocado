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
package org.bdgenomics.avocado.cli

import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.AvocadoFunSuite

class MergeDiscoveredSuite extends AvocadoFunSuite {

  sparkTest("merge variants discovered from two samples") {
    val path = resourceUrl("NA12878.1_907170_A2C.vcf")
      .toString
      .replace("A2C", "*")
    val loc = tmpLocation()
    MergeDiscovered(Array(path, loc)).run(sc)
    assert(sc.loadVariants(loc).dataset.count === 3)
  }
}
