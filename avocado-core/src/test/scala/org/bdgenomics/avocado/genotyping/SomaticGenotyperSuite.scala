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
package org.bdgenomics.avocado.genotyping

import org.bdgenomics.adam.models.{
  ReferencePosition,
  SequenceDictionary,
  SequenceRecord,
  VariantContext
}
import org.bdgenomics.avocado.models.{ Observation, AlleleObservation }
import org.bdgenomics.formats.avro.Contig
import org.scalatest.FunSuite
import scala.collection.JavaConversions._
import scala.math.{ abs, sqrt }

class SomaticGenotyperSuite extends FunSuite {
  val sg = new SomaticGenotyper(SequenceDictionary(SequenceRecord("ctg", 1000L)),
    "normalSample",
    "somaticSample")

  test("don't call a variant with a hom-ref normal and hom-ref somatic") {
    val observed = Iterable(
      new Observation(ReferencePosition("ctg", 0L), "C"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        40,
        true,
        "somaticSample"))

    val call = sg.genotypeSite((ReferencePosition("ctg", 0L), observed))

    assert(call.size === 0)
  }

  test("call variant for hom-ref normal, \"hom-alt\" somatic") {
    val observed = Iterable(
      new Observation(ReferencePosition("ctg", 0L), "C"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        30,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        40,
        40,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        40,
        true,
        "somaticSample"))

    val call = sg.genotypeSite((ReferencePosition("ctg", 0L), observed))

    assert(call.size === 1)
    assert(call.head.variant.getAlternateAllele === "A")
    assert(call.head.genotypes.size === 1)
    assert(call.head.genotypes.head.getSampleId === "somaticSample")
    assert(call.head.genotypes.head.getGenotypeQuality === 105)
    assert(call.head.genotypes.head.getReadDepth === 3)
    assert(call.head.genotypes.head.getAlternateReadDepth === 3)
  }

  test("don't call a variant for matching het normal and somatic") {
    val observed = Iterable(
      new Observation(ReferencePosition("ctg", 0L), "C"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        40,
        40,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        40,
        true,
        "somaticSample"))

    val call = sg.genotypeSite((ReferencePosition("ctg", 0L), observed))

    assert(call.size === 0)
  }

  test("call a variant for het normal and tri-allelic somatic") {
    val observed = Iterable(
      new Observation(ReferencePosition("ctg", 0L), "C"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        40,
        40,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "T",
        30,
        40,
        true,
        "somaticSample"))

    val call = sg.genotypeSite((ReferencePosition("ctg", 0L), observed))

    assert(call.size === 1)
    assert(call.head.variant.getAlternateAllele === "T")
    assert(call.head.genotypes.size === 1)
    assert(call.head.genotypes.head.getSampleId === "somaticSample")
    assert(call.head.genotypes.head.getGenotypeQuality === 35)
    assert(call.head.genotypes.head.getReadDepth === 3)
    assert(call.head.genotypes.head.getAlternateReadDepth === 1)
  }

  test("call two variants for hom-ref normal and tri-allelic somatic") {
    val observed = Iterable(
      new Observation(ReferencePosition("ctg", 0L), "C"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        40,
        true,
        "normalSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        40,
        40,
        true,
        "somaticSample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "T",
        30,
        40,
        true,
        "somaticSample"))

    val call = sg.genotypeSite((ReferencePosition("ctg", 0L), observed))

    assert(call.size === 2)
    assert(call.filter(_.variant.getAlternateAllele == "A").size === 1)
    assert(call.filter(_.variant.getAlternateAllele == "T").size === 1)
    assert(call.filter(_.variant.getAlternateAllele == "A").head.variant.getAlternateAllele === "A")
    assert(call.filter(_.variant.getAlternateAllele == "A").head.genotypes.size === 1)
    assert(call.filter(_.variant.getAlternateAllele == "A").head.genotypes.head.getSampleId === "somaticSample")
    assert(call.filter(_.variant.getAlternateAllele == "A").head.genotypes.head.getGenotypeQuality === 40)
    assert(call.filter(_.variant.getAlternateAllele == "A").head.genotypes.head.getReadDepth === 3)
    assert(call.filter(_.variant.getAlternateAllele == "A").head.genotypes.head.getAlternateReadDepth === 1)
    assert(call.filter(_.variant.getAlternateAllele == "T").head.variant.getAlternateAllele === "T")
    assert(call.filter(_.variant.getAlternateAllele == "T").head.genotypes.size === 1)
    assert(call.filter(_.variant.getAlternateAllele == "T").head.genotypes.head.getSampleId === "somaticSample")
    assert(call.filter(_.variant.getAlternateAllele == "T").head.genotypes.head.getGenotypeQuality === 35)
    assert(call.filter(_.variant.getAlternateAllele == "T").head.genotypes.head.getReadDepth === 3)
    assert(call.filter(_.variant.getAlternateAllele == "T").head.genotypes.head.getAlternateReadDepth === 1)
  }
}
