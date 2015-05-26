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
import org.bdgenomics.avocado.algorithms.math.MathTestUtils
import org.bdgenomics.avocado.discovery.ReadExplorer
import org.bdgenomics.avocado.models.{ Observation, AlleleObservation }
import org.bdgenomics.formats.avro.{ AlignmentRecord, Contig }
import org.scalatest.FunSuite
import scala.collection.JavaConversions._
import scala.math.{ log, sqrt }

class BiallelicGenotyperSuite extends FunSuite {
  val ba = new BiallelicGenotyper(SequenceDictionary(SequenceRecord("ctg", 1000L)),
    Map(("ctg", 1000L)))

  test("score genotype for single sample, all bases ref") {
    val observed = Iterable(
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        Some(30),
        true,
        true,
        1,
        "mySample",
        1L),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        Some(40),
        true,
        true,
        1,
        "mySample",
        2L),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        Some(40),
        true,
        true,
        1,
        "mySample",
        3L))

    val expected = Array(8.0 * (0.001 *
      0.0001 *
      sqrt(0.001 * 0.0001)) / 8.0,
      0.125,
      8.0 * (0.999 *
        0.9999 *
        (1.0 - sqrt(0.001 * 0.0001))) / 8.0).map(log(_))

    val scored = ba.scoreGenotypeLikelihoods("C", "A", observed)

    MathTestUtils.assertAlmostEqual(expected(0), scored._2(0))
    MathTestUtils.assertAlmostEqual(expected(1), scored._2(1))
    MathTestUtils.assertAlmostEqual(expected(2), scored._2(2))
    assert(scored._2.max == scored._2(2))
  }

  test("score genotype for single sample, mix of ref/non-ref bases") {
    val observed = Iterable(
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        Some(30),
        true,
        true,
        1,
        "mySample",
        1L),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        Some(40),
        true,
        true,
        1,
        "mySample",
        2L),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        Some(40),
        true,
        true,
        1,
        "mySample",
        3L))

    val expected = List(8.0 * (0.001 *
      0.0001 *
      (1.0 - sqrt(0.001 * 0.0001))) / 8.0,
      0.125,
      8.0 * (0.999 *
        0.9999 *
        sqrt(0.001 * 0.0001)) / 8.0).map(log(_))

    val scored = ba.scoreGenotypeLikelihoods("C", "A", observed)

    MathTestUtils.assertAlmostEqual(expected(0), scored._2(0))
    MathTestUtils.assertAlmostEqual(expected(1), scored._2(1))
    MathTestUtils.assertAlmostEqual(expected(2), scored._2(2))
    assert(scored._2.max == scored._2(1))
  }

  test("score genotype for single sample, all bases non-ref") {
    val observed = Iterable(
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        Some(30),
        true,
        true,
        1,
        "mySample",
        1L),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        40,
        Some(40),
        true,
        true,
        1,
        "mySample",
        2L),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        Some(40),
        true,
        true,
        1,
        "mySample",
        3L))

    val expected = List(8.0 * (0.999 *
      0.9999 *
      (1.0 - sqrt(0.001 * 0.0001))) / 8.0,
      0.125,
      8.0 * (0.001 *
        0.0001 *
        sqrt(0.001 * 0.0001)) / 8.0).map(log(_))

    val scored = ba.scoreGenotypeLikelihoods("C", "A", observed)

    MathTestUtils.assertAlmostEqual(expected(0), scored._2(0))
    MathTestUtils.assertAlmostEqual(expected(1), scored._2(1))
    MathTestUtils.assertAlmostEqual(expected(2), scored._2(2))
    assert(scored._2.max == scored._2(0))
  }

  /**
   * test("emit genotypes on iterator") {
   * val re = new ReadExplorer(null)
   *
   * val read = AlignmentRecord.newBuilder()
   * .setStart(10L)
   * .setEnd(15L)
   * .setContig(Contig.newBuilder()
   * .setContigName("ctg")
   * .build())
   * .setMapq(40)
   * .setSequence("ACTGA")
   * .setQual(":::::")
   * .setCigar("1M1I2M1D1M")
   * .setRecordGroupSample("sample1")
   * .build()
   *
   * val observations = re.readToObservations((read, 0L)).toSeq ++ Seq(
   * new Observation(ReferencePosition("ctg", 10L), "A", 1),
   * new Observation(ReferencePosition("ctg", 11L), "T", 1),
   * new Observation(ReferencePosition("ctg", 12L), "G", 1),
   * new Observation(ReferencePosition("ctg", 13L), "T", 1),
   * new Observation(ReferencePosition("ctg", 14L), "A", 1))
   *
   * val gts = ba.genotypeIterator(observations.sortBy(_.pos).map(v => (v.pos, v)).toIterator).toSeq
   *
   * assert(gts.length === 4)
   * gts.filter(_.position.pos != 12L).map(_.variant).foreach(v => {
   * assert(v.getReferenceAllele.length === 1)
   * })
   * gts.filter(_.position.pos == 12L).map(_.variant).foreach(v => {
   * assert(v.getReferenceAllele.length === 2)
   * })
   * gts.filter(_.position.pos != 10L).map(_.variant).foreach(v => {
   * assert(v.getAlternateAllele.length === 1)
   * })
   * gts.filter(_.position.pos == 10L).map(_.variant).foreach(v => {
   * assert(v.getAlternateAllele.length === 2)
   * })
   * }
   */
}
