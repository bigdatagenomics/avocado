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

package edu.berkeley.cs.amplab.avocado.filters.pileup

import spark.{RDD,SparkContext}
import edu.berkeley.cs.amplab.adam.util.{Pileup,PileupTraversable}

/**
 * Trait for filtering pileups. 
 */
trait PileupFilter {

  /**
   * Method signature for filter operation.
   *
   * @param[in] pileups An RDD containing reference oriented stacks of nucleotides.
   * @return An RDD containing lists of pileups.
   */
  def filter (pileups: RDD [(void, Pileup)]): RDD [(Any, List[Pileup])] 
}
