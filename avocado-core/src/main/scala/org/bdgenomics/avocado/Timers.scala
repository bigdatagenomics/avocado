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
 * Contains [[Timers]] that are used to instrument ADAM.
 */
private[avocado] object Timers extends Metrics {

  // File Loading
  val LoadReads = timer("Load Reads")
  val LoadContigs = timer("Load Reference Contigs")

  // Pipeline Stages
  val PreprocessReads = timer("Preprocessing Reads")
  val CallVariants = timer("Call Variants")
  val DiscoverObservations = timer("Discover Observations")
  val GenotypeObservations = timer("Genotype Observations")
  val PostprocessVariants = timer("Postprocess Variants")

  // Biallelic Genotyping Model
  val GenotypingSite = timer("Scoring Genotypes at a Site")
  val ScoringLikelihoods = timer("Scoring Likelihoods")
  val EmittingCall = timer("Assembling Genotype Call and Statistics")

  // Simple Read Explorer
  val ExploringReads = timer("Extracting Read Observations")
  val ExploringRead = timer("Observing a Read")

  // Reassembly Explorer
  val JoiningReads = timer("Joining Reads and References Together")
  val ProcessingRegions = timer("Processing Regions")
  val RegionDiscovery = timer("Observing Region")
  val CheckActivity = timer("Calculating Region Activity Statistics")
  val ReassemblingRegion = timer("Reassembling Region")
  val BuildingGraph = timer("Building Indexed de Bruijn Graph")
  val ObservingGraph = timer("Crawling Indexed de Bruijn Graph")
  val InactiveReads = timer("Extracting Read Observations from Inactive Regions")

  // Statistics
  val ComputingCoverage = timer("Computing Coverage")
  val ReferenceLengths = timer("Computing Contig Lengths")
  val CollectingReference = timer("Collecting Reference Sequence")
  val ExtractingSequenceDictionary = timer("Extracting Sequence Dictionary")
  val CollectingSamples = timer("Collecting Sample Names")
  val ExploringReference = timer("Extracting Reference Observations")

  // File Saving
  val SaveVariants = timer("Save Variants")
}
