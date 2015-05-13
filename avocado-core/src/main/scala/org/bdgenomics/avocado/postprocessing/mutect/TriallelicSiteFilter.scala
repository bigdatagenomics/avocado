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
 * Based on the following text:
 *   Reject false positives caused by calling triallelic sites where the normal
 *   sample is heterozygous with alleles A and B, and MuTect is considering an
 *   alternate allele C. Although this is biologically possible, and remains an
 *   area for future improvement in the detection of mutations, calling at these
 *   sites generates many false positives, and therefore they are currently
 *   filtered out by default. However, it may be desirable to review mutations
 *   failing only this filter for biological relevance and orthogonal validation,
 *   and to study the underlying reasons for these false positives.
 */
class TriallelicSiteFilter extends MutectPostprocessor with Serializable {

  def isTriAllelic(vc: VariantContext, normalReads: Iterable[AlignmentRecord]): Boolean =
    ???

  override def filter(variants: RDD[VariantContext],
                      tumorReads: RDD[Classified[AlignmentRecord]],
                      normalReads: RDD[Classified[AlignmentRecord]]): RDD[VariantContext] = {
    ???
  }
}
