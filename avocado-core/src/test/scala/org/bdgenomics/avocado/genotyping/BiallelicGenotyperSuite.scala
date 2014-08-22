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

import org.bdgenomics.adam.models.{ ReferencePosition, VariantContext }
import org.bdgenomics.avocado.models.{ Observation, AlleleObservation }
import org.scalatest.FunSuite
import scala.collection.JavaConversions._
import scala.math.{ abs, sqrt }

class BiallelicGenotyperSuite extends FunSuite {
  val ba = new BiallelicGenotyper()
  val floatingPointingThreshold = 1e-6

  def assertAlmostEqual(a: Double, b: Double, epsilon: Double = floatingPointingThreshold) {
    if (!(a * 0.99 < b && a * 1.01 > b) &&
      !(abs(a - b) < epsilon)) {
      throw new AssertionError(a + " != " + b)
    }
  }

  test("score genotype for single sample, all bases ref") {
    val observed = Iterable(
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "mySample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "mySample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        40,
        true,
        "mySample"))

    val expected = Array(8.0 * (0.001 *
      0.0001 *
      sqrt(0.001 * 0.0001)) / 8.0,
      0.125,
      8.0 * (0.999 *
        0.9999 *
        (1.0 - sqrt(0.001 * 0.0001))) / 8.0)

    val scored = ba.scoreGenotypeLikelihoods("C", "A", observed)

    assertAlmostEqual(expected(0), scored._2(0))
    assertAlmostEqual(expected(1), scored._2(1))
    assertAlmostEqual(expected(2), scored._2(2))
    assert(scored._2.max == scored._2(2))
  }

  test("score genotype for single sample, mix of ref/non-ref bases") {
    val observed = Iterable(
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        30,
        30,
        true,
        "mySample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "C",
        40,
        40,
        true,
        "mySample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        40,
        true,
        "mySample"))

    val expected = List(8.0 * (0.001 *
      0.0001 *
      (1.0 - sqrt(0.001 * 0.0001))) / 8.0,
      0.125,
      8.0 * (0.999 *
        0.9999 *
        sqrt(0.001 * 0.0001)) / 8.0)

    val scored = ba.scoreGenotypeLikelihoods("C", "A", observed)

    assertAlmostEqual(expected(0), scored._2(0))
    assertAlmostEqual(expected(1), scored._2(1))
    assertAlmostEqual(expected(2), scored._2(2))
    assert(scored._2.max == scored._2(1))
  }

  test("score genotype for single sample, all bases non-ref") {
    val observed = Iterable(
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        30,
        true,
        "mySample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        40,
        40,
        true,
        "mySample"),
      AlleleObservation(ReferencePosition("ctg", 0L),
        1,
        "A",
        30,
        40,
        true,
        "mySample"))

    val expected = List(8.0 * (0.999 *
      0.9999 *
      (1.0 - sqrt(0.001 * 0.0001))) / 8.0,
      0.125,
      8.0 * (0.001 *
        0.0001 *
        sqrt(0.001 * 0.0001)) / 8.0)

    val scored = ba.scoreGenotypeLikelihoods("C", "A", observed)

    assertAlmostEqual(expected(0), scored._2(0))
    assertAlmostEqual(expected(1), scored._2(1))
    assertAlmostEqual(expected(2), scored._2(2))
    assert(scored._2.max == scored._2(0))
  }
}
