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
 *   Eliminate false positives in the tumor data by looking at the control data
 *   (typically from the matched normal sample) for evidence of the alternate
 *   allele beyond what is expected from random sequencing error. A candidate
 *   is rejected if, in the control data, there are
 *     (i) ≥ 2 observations of the alternate allele or they represent ≥ 3% of
 *         the reads; and
 *     (ii) their sum of quality scores is > 20.
 */
class ObservedInControlFilter extends MutectPostprocessor {
  override def filter(variants: RDD[VariantContext],
                      tumorReads: RDD[Classified[AlignmentRecord]],
                      normalReads: RDD[Classified[AlignmentRecord]]): RDD[VariantContext] =
    ???
}
