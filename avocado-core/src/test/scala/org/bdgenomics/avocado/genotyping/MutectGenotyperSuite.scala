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
  val allC = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = "C",
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID * 10,
      100,
      None, None,
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMuts = for (readID <- 1 to 20) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID * 10 % 100,
      100,
      None, Some(-50),
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val allMuts = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = "A", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID * 10,
      100,
      None, None,
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsEndClusteredBeginning = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID,
      100,
      None, None,
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsSoftClipped = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID,
      100,
      None, None,
      20, 20, 100, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsNoisyReads = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID,
      100,
      None, None,
      0, 0, 100, Some(101), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsNearInsertion = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID,
      100,
      Some(3), None,
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsNearDeletion = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID,
      100,
      None, Some(3),
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsMateRescue = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID,
      100,
      None, Some(3),
      0, 0, 1, Some(0), true,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsEndClusteredEnd = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = 100 - readID,
      100,
      None, None,
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsStrandBiasNeg = for (readID <- 1 to 50) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 25) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID <= 25,
      firstOfPair = true,
      offsetInAlignment = readID * 10 % 100,
      100,
      None, Some(-50),
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val someMutsStrandBiasPos = for (readID <- 1 to 50) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID <= 25) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID > 25,
      firstOfPair = true,
      offsetInAlignment = readID * 10 % 100,
      100,
      None, Some(-50),
      0, 0, 1, Some(0), false,
      sample = "tumor",
      readId = 1L)
  }

  val normalAllC = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = "C",
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID * 10,
      100,
      None, None,
      0, 0, 1, Some(0), false,
      sample = "normal",
      readId = 1L)
  }
  val normalHet = for (readID <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (readID % 2 == 0) "C" else "A",
      phred = 30 + (readID - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = readID % 2 == 0,
      firstOfPair = true,
      offsetInAlignment = readID * 10,
      100,
      None, None,
      0, 0, 1, Some(0), false,
      sample = "normal",
      readId = 1L)
  }

  val noMutClean = allC ++ normalAllC
  val noMutHet = allC ++ normalHet
  val mutClean = someMuts ++ normalAllC
  val mutHet = someMuts ++ normalHet
  val mutEndClusterBeginning = someMutsEndClusteredBeginning ++ normalAllC
  val mutEndClusterEnd = someMutsEndClusteredEnd ++ normalAllC
  val mutHeavilyClipped = someMutsSoftClipped ++ normalAllC
  val mutNearInsertion = someMutsNearInsertion ++ normalAllC
  val mutNearDeletion = someMutsNearDeletion ++ normalAllC
  val mutNoisyReads = someMutsNoisyReads ++ normalAllC
  val mutMateRescue = someMutsMateRescue ++ normalAllC
  val mutStrandBiasNeg = someMutsStrandBiasNeg ++ normalAllC
  val mutStrandBiasPos = someMutsStrandBiasPos ++ normalAllC
  val mutAllMutants = allMuts ++ normalAllC

  test("Sites with no mutations should not have any variants returned") {
    val resultClean = mt.genotypeSite(pos, ref, noMutClean)
    assert(resultClean.isEmpty, "Result should be empty")
    val resultHet = mt.genotypeSite(pos, ref, noMutHet)
    assert(resultHet.isEmpty, "Although normal is het, tumor has no mutants and there should be no result")
    val resultEmptyTumor = mt.genotypeSite(pos, ref, normalHet)
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
    val resultMissingNormal = mt.genotypeSite(pos, ref, someMuts)
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

  test("Test strand bias filter") {
    val resultNeg = mt.genotypeSite(pos, ref, mutStrandBiasNeg)
    val resultPos = mt.genotypeSite(pos, ref, mutStrandBiasPos)
    assert(resultNeg.isEmpty, "Fails negative strand bias filter")
    assert(resultPos.isEmpty, "Fails positive strand bias filter")
  }

  test("Works even when all alleles are mutant") {
    val result = mt.genotypeSite(pos, ref, mutAllMutants)
    assert(result.isDefined)
  }

}
