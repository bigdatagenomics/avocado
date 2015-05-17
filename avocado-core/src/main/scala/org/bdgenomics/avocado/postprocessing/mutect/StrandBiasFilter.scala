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
 * Text from the Mutect paper:
 *
 *   Reject false positives caused by context-specific sequencing errors where the vast
 *   majority of the alternate alleles are observed in a single direction of reads. We
 *   perform this test by stratifying the reads by direction and then applying the core
 *   detection statistic on the two data sets. We also calculate the sensitivity to
 *   have passed the threshold given the data (Online Methods). Candidates are rejected
 *   when the strand-specific LOD is < 2.0 in directions where the sensitivity to have
 *   passed that threshold is ≥ 90%.
 */
class StrandBiasFilter extends MutectPostprocessor {
  override def filter(variants: RDD[VariantContext],
                      tumorReads: RDD[Classified[AlignmentRecord]],
                      normalReads: RDD[Classified[AlignmentRecord]]): RDD[VariantContext] = {
    ???
  }
}
