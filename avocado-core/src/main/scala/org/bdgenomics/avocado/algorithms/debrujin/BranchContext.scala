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
package org.bdgenomics.avocado.algorithms.debrujin

import org.bdgenomics.adam.models.ReferencePosition

private[debrujin] trait BranchContext {
  val kmer: Kmer
}

object Allele {
  def apply(kmer: Kmer,
            branchPoint: ReferencePosition): Allele = {
    Allele(kmer,
      "",
      branchPoint,
      List(),
      kmer.readId.toSet)
  }
}

private[debrujin] case class Allele(kmer: Kmer,
                                    allele: String,
                                    branchPoint: ReferencePosition,
                                    pending: List[Kmer],
                                    activeReads: Set[Long]) extends BranchContext {
  assert(!kmer.isReference, "Cannot create allele for mapped k-mer. Must be closed allele.")
}

private[debrujin] case class ClosedAllele(kmer: Kmer,
                                          allele: String,
                                          branchPoint: ReferencePosition,
                                          pending: List[Kmer],
                                          activeReads: Set[Long]) extends BranchContext {
  assert(kmer.isReference, "K-mer must be mapped to close allele.")
}

private[debrujin] case class Reference(kmer: Kmer) extends BranchContext {
  assert(kmer.isReference, "K-mer must be mapped to a reference position.")
}
