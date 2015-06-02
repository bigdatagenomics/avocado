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

package org.bdgenomics.avocado.genotyping

import org.bdgenomics.adam.models.ReferencePosition
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.scalatest.FunSuite

class MutectGenotyperSuite extends FunSuite {
  val pos = ReferencePosition("ctg", 101l)
  val mt = new MutectGenotyper("normal", "tumor")

  val ref = new Observation(ReferencePosition("ctg", 101L), "C")

  //Demo data
  val all_c = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = "C",
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id * 10,
      100,
      None, None,
      0, 0, 1, 0, false,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id * 10,
      100,
      None, Some(-50),
      0, 0, 1, 0, false,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts_end_clustered_beginning = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id,
      100,
      None, None,
      0, 0, 1, 0, false,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts_soft_clipped = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id,
      100,
      None, None,
      20, 20, 100, 0, false,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts_noisy_reads = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id,
      100,
      None, None,
      0, 0, 100, 101, false,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts_near_insertion = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id,
      100,
      Some(3), None,
      0, 0, 1, 0, false,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts_near_deletion = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id,
      100,
      None, Some(3),
      0, 0, 1, 0, false,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts_mate_rescue = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id,
      100,
      None, Some(3),
      0, 0, 1, 0, true,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts_end_clustered_end = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = 100 - read_id,
      100,
      None, None,
      0, 0, 1, 0, false,
      sample = "tumor",
      readId = 1L)
  }

  val normal_all_c = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = "C",
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id * 10,
      100,
      None, None,
      0, 0, 1, 0, false,
      sample = "normal",
      readId = 1L)
  }
  val normal_het = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id % 2 == 0) "C" else "A",
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = read_id % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = read_id * 10,
      100,
      None, None,
      0, 0, 1, 0, false,
      sample = "normal",
      readId = 1L)
  }

  val noMutClean = all_c ++ normal_all_c
  val noMutHet = all_c ++ normal_het
  val mutClean = some_muts ++ normal_all_c
  val mutHet = some_muts ++ normal_het
  val mutEndClusterBeginning = some_muts_end_clustered_beginning ++ normal_het
  val mutEndClusterEnd = some_muts_end_clustered_end ++ normal_het
  val mutHeavilyClipped = some_muts_soft_clipped ++ normal_het
  val mutNearInsertion = some_muts_near_insertion ++ normal_het
  val mutNearDeletion = some_muts_near_deletion ++ normal_het
  val mutNoisyReads = some_muts_noisy_reads ++ normal_het
  val mutMateRescue = some_muts_mate_rescue ++ normal_het

  test("Sites with no mutations should not have any variants returned") {
    val resultClean = mt.genotypeSite(pos, ref, noMutClean)
    assert(resultClean.isEmpty, "Result should be empty")
    val resultHet = mt.genotypeSite(pos, ref, noMutHet)
    assert(resultHet.isEmpty, "Although normal is het, tumor has no mutants and there should be no result")
    val resultEmptyTumor = mt.genotypeSite(pos, ref, normal_het)
    assert(resultEmptyTumor.isEmpty, "A position with no tumor reads should not be classified somatic.")
  }

  test("Sites with mutations should return a variant") {
    val result = mt.genotypeSite(pos, ref, mutClean)
    assert(result.isDefined, "Result should not be empty")
    val variant = result.get.variant.variant
    assert(variant.getAlternateAllele === "A", "Alternate allele should be a")
    assert(variant.getReferenceAllele === "C", "Reference allele should be c")
    assert(variant.getContig.getContigName === "ctg", "Contig should be 'ctg'")
    assert(variant.getStart === 101l, "Variant should be at position 101")
    assert(variant.getEnd === 102l, "Variant should be at position 101")
  }

  test("Sites with tumor mutations that are heterozygous in normal should not return a variant") {
    val result = mt.genotypeSite(pos, ref, mutHet)
    assert(result.isEmpty, "Site should not be classified somatic, strong het in normal")
    val resultMissingNormal = mt.genotypeSite(pos, ref, some_muts)
    assert(resultMissingNormal.isEmpty, "When there is no normal coverage, the site should not be classified somatic")
  }

  test("Mutations that cluster near beginning/end of alignment should be discarded") {
    val beginning = mt.genotypeSite(pos, ref, mutEndClusterBeginning)
    val end = mt.genotypeSite(pos, ref, mutEndClusterEnd)
    assert(beginning.isEmpty, "Mutation clustering at beginning of the alignment should fail mutect")
    assert(end.isEmpty, "Mutation clustering at end of the alignment should fail mutect")

  }

  test("Test noisy read filter") {
    val result = mt.genotypeSite(pos, ref, mutNoisyReads)
    assert(result.isEmpty)
  }

  test("Test heavily clipped filter") {
    val result = mt.genotypeSite(pos, ref, mutHeavilyClipped)
    assert(result.isEmpty)
  }

  test("Test near insertion filter") {
    val result = mt.genotypeSite(pos, ref, mutNearInsertion)
    assert(result.isEmpty)
  }

  test("Test near deletion filter") {
    val result = mt.genotypeSite(pos, ref, mutNearDeletion)
    assert(result.isEmpty)
  }

  test("Test mate rescue filter") {
    val result = mt.genotypeSite(pos, ref, mutMateRescue)
    assert(result.isEmpty)
  }

}
