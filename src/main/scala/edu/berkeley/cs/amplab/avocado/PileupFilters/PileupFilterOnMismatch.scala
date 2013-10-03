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
import edu.berkeley.cs.amplab.avocado.filters.pileup.PileupFilter

/**
 * Class implementing a filter that only returns Pileups that contain
 * a mismatch on any bases in the pileup.
 */
class PileupFilterOnMismatch extends PileupFilter {

  /**
   * Filter to only return pileups with mismatches.
   *
   * @param[in] pileups An RDD containing reference oriented stacks of nucleotides.
   * @return An RDD containing only pileups that contain at least one mismatch.
   */
  override def filter (pileups: RDD [(void, Pileup)]): RDD [(void, List[Pileup])] = {
    return pileups.filter (_.mismatches != List.empty)
                  .map (_ => (void, List (_)))
  }
}
