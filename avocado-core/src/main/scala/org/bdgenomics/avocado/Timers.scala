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
package org.bdgenomics.avocado

import org.bdgenomics.utils.instrumentation.Metrics

/**
 * Contains [[Timers]] that are used to instrument Avocado.
 */
private[avocado] object Timers extends Metrics {

  // org.bdgenomics.avocado.genotyping.BiallelicGenotyper
  val DiscoverAndCall = timer("Discover variants and call genotypes")
  val CallGenotypes = timer("Call genotypes")
  val JoinReadsAndVariants = timer("Joining reads against variants")
  val ObserveReads = timer("Observing variants in reads")
  val ObserveRead = timer("Observing a read")
  val IntersectVariants = timer("Intersecting observations and variants")
  val ProcessIntersections = timer("Processing intersections")
  val ProcessException = timer("Processing exceptions")
  val EmitGenotypes = timer("Emit observations as genotype calls")

  // org.bdgenomics.avocado.genotyping.DiscoverVariants
  val DiscoveringVariants = timer("Discovering variants in reads")

  // org.bdgenomics.avocado.models.ObservationOperator
  val ExtractingReference = timer("Extracting reference sequence")
  val ExtractingAlignment = timer("Parsing alignment (CIGAR/MD tag)")

  // org.bdgenomics.avocado.realigner.Aligner
  val AligningSequences = timer("Aliging read against reference")

  // org.bdgenomics.avocado.realigner.Realigner
  val ProcessingReadForRealignment = timer("Processing a single read for realignment")
  val RealigningRead = timer("Realigning candidate read")
  val FinalizingRealignment = timer("Finalizing realignment")

  // org.bdgenomics.avocado.realigner.RealignmentBlock
  val ExtractingRealignmentBlocks = timer("Checking for non-canonical alignment blocks")

  // org.bdgenomics.avocado.util.Downsampler
  val DownsampleReads = timer("Downsampling reads across an RDD")
  val DownsamplePartition = timer("Downsampling a single partition")

  // org.bdgenomics.avocado.util.GapFilter
  val FilteringByGaps = timer("Filtering out reads aligned to gaps")

  // org.bdgenomics.avocado.util.TreeRegionJoin
  val TreeJoin = timer("Running broadcast join with interval tree")
  val BuildingTrees = timer("Building interval tree")
  val SortingRightSide = timer("Sorting right side of join")
  val GrowingTrees = timer("Growing forest of trees")
  val RunningMapSideJoin = timer("Running map-side join")
}
