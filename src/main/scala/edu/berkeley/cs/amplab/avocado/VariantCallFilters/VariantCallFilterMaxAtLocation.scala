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

package edu.berkeley.cs.amplab.avocado.filters.variants

import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.{ADAMVariant}

/**
 * Simple filter. Checks to make sure that no more variants show up at a
 * single location than the ploidy of the organism.
 *
 * For now, if more than two variants show up at a single location, we throw
 * those variants out. This will be changed later.
 */
class VariantCallFilterMaxAtLocation extends VariantCallFilter {

  val filterName = "MaxAtLocation"

  // get ploidy - if no ploidy, assume human diploid
  val ploidy = 2
    
  /**
   * Simple filter to check for no more than two variants at a single location.
   *
   * @param[in] variants An RDD containing variants.
   * @return An RDD containing variants.
   */
  override def filter (variants: RDD [ADAMVariant]): RDD [ADAMVariant] = {
    variants.groupBy (_.getStartPosition).filter (kv => kv._2.length <= ploidy)
      .flatMap (kv => kv._2.toList)
  }

}
