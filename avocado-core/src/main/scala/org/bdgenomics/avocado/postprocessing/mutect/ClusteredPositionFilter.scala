/*
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.avocado.postprocessing.mutect

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.formats.avro.AlignmentRecord

/**
 *   Reject false positives caused by mis-alignments hallmarked by the alternate
 *   alleles being clustered at a consistent distance from the start or end of
 *   the read alignment. We calculate the median and median absolute deviation
 *   of the distance from both the start and end of the read and reject sites
 *   that have a median ≤ 10 (near the start/end of the alignment) and a median
 *   absolute deviation ≤ 3 (clustered).
 */
class ClusteredPositionFilter extends MutectPostprocessor {
  override def filter(variants: RDD[VariantContext],
                      tumorReads: RDD[Classified[AlignmentRecord]],
                      normalReads: RDD[Classified[AlignmentRecord]]): RDD[VariantContext] =
    ???
}
