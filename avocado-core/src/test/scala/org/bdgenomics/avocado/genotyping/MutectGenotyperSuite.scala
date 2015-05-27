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
      offsetInRead = 1,
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
      offsetInRead = 1,
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
      offsetInRead = 1,
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
      offsetInRead = 1,
      sample = "normal",
      readId = 1L)
  }

  val noMutClean = all_c ++ normal_all_c
  val noMutHet = all_c ++ normal_het
  val mutClean = some_muts ++ normal_all_c
  val mutHet = some_muts ++ normal_het

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
}
