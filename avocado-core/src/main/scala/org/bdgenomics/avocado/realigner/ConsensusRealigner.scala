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
package org.bdgenomics.avocado.realigner

import org.bdgenomics.adam.algorithms.consensus.ConsensusGenerator
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variant.VariantRDD
import org.bdgenomics.avocado.genotyping.DiscoverVariants
import org.bdgenomics.formats.avro.Variant

object ConsensusRealigner {

  /**
   * Realigns a set of reads against the reference genome.
   *
   * @param reads Reads to realign.
   * @param kmerLength The length k of the k-mers.
   * @return Returns the realigned reads.
   */
  def realign(reads: AlignmentRecordRDD,
              kmerLength: Int): AlignmentRecordRDD = {

    // use the realigner to make a first pass over the reads
    val realignedReads = Realigner.realign(reads, kmerLength)

    // from here, we'll discover any potential variants
    val variants = filterIndels(DiscoverVariants(realignedReads))

    // we'll pass these discovered variants to ADAM's indel realigner
    reads.realignIndels(ConsensusGenerator.fromKnownIndels(variants))
  }

  /**
   * @param v The variant to filter.
   * @return Returns false if the variant is a SNV/MNV.
   */
  private[realigner] def discardSnvs(v: Variant): Boolean = {
    v.getReferenceAllele.length != v.getAlternateAllele.length
  }

  /**
   * @param rdd The RDD of variants to filter.
   * @return Returns a new RDD of variants where the SNVs & MNVs have been
   *   removed and only INDEL variants remain.
   */
  private[realigner] def filterIndels(rdd: VariantRDD): VariantRDD = {
    rdd.transform(r => r.filter(discardSnvs))
  }
}
