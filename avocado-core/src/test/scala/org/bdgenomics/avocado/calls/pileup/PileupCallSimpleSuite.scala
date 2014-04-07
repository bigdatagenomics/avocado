/*
 * Copyright (c) 2013. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.avocado.calls.pileup

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.avro.{ ADAMPileup, Base }
import org.bdgenomics.adam.models.{ ADAMRod, ADAMVariantContext }
import scala.collection.JavaConversions._
import org.scalatest.FunSuite
import scala.math.abs

class PileupCallSimpleSuite extends FunSuite {

  val floatingPointingThreshold = 1e-6

  def assertAlmostEqual(a: Double, b: Double, epsilon: Double = floatingPointingThreshold) {
    assert((a * 0.99 < b && a * 1.01 > b) ||
      abs(a - b) < epsilon)
  }

  test("Collect the max non-ref base for homozygous SNP") {
    val pl = List(ADAMPileup.newBuilder()
      .setReadBase(Base.A)
      .setReferenceBase(Base.C)
      .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.A)
        .setReferenceBase(Base.C)
        .build())

    val call = new PileupCallSimpleSNP(2)

    val mb = call.getMaxNonRefBase(pl)

    assert(mb.isDefined)
    assert(mb.get === Base.A)
  }

  test("Collect the max non-ref base for hetereozygous SNP") {
    val pl = List(ADAMPileup.newBuilder()
      .setReadBase(Base.C)
      .setReferenceBase(Base.C)
      .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.A)
        .setReferenceBase(Base.C)
        .build())

    val call = new PileupCallSimpleSNP(2)

    val mb = call.getMaxNonRefBase(pl)

    assert(mb.isDefined)
    assert(mb.get === Base.A)
  }

  test("Collect the max non-ref base for homozygous SNP that is not biallelic") {
    val pl = List(ADAMPileup.newBuilder()
      .setReadBase(Base.A)
      .setReferenceBase(Base.C)
      .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.A)
        .setReferenceBase(Base.C)
        .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.T)
        .setReferenceBase(Base.C)
        .build())

    val call = new PileupCallSimpleSNP(2)

    val mb = call.getMaxNonRefBase(pl)

    assert(mb.isDefined)
    assert(mb.get === Base.A)
  }

  test("Do not collect the max non-ref base for homozygous ref.") {
    val pl = List(ADAMPileup.newBuilder()
      .setReadBase(Base.C)
      .setReferenceBase(Base.C)
      .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.C)
        .setReferenceBase(Base.C)
        .build())

    val call = new PileupCallSimpleSNP(2)

    val mb = call.getMaxNonRefBase(pl)

    assert(mb.isEmpty)
  }

  test("Default compensation method doesn't perform any compensation") {
    val call = new PileupCallSimpleSNP(2)

    val l = List(0.0, 0.0, 0.0)

    assert(l === call.compensate(l, List[ADAMPileup]()))
  }

  test("score genotype for single sample, all bases ref") {
    val call = new PileupCallSimpleSNP(2)

    val pl = List(ADAMPileup.newBuilder()
      .setReadBase(Base.C)
      .setReferenceBase(Base.C)
      .setMapQuality(30)
      .setSangerQuality(30)
      .setCountAtPosition(1)
      .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.C)
        .setReferenceBase(Base.C)
        .setMapQuality(40)
        .setSangerQuality(40)
        .setCountAtPosition(1)
        .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.C)
        .setReferenceBase(Base.C)
        .setMapQuality(40)
        .setSangerQuality(30)
        .setCountAtPosition(1)
        .build())

    val expected = List((8.0 * (0.999 * 0.999 * 0.9999 * 0.9999 * 0.999 * 0.9999)) / 8.0,
      ((0.999 * 0.999 * 0.9999 * 0.9999 * 0.999 * 0.9999) +
        (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999)) / 8.0,
      8.0 * (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999) / 8.0).reverse

    val scored = call.scoreGenotypeLikelihoods(pl)

    for (i <- 0 to 2) {
      assertAlmostEqual(expected(i), scored(i))
    }
  }

  test("score genotype for single sample, mix of ref/non-ref bases") {
    val call = new PileupCallSimpleSNP(2)

    val pl = List(ADAMPileup.newBuilder()
      .setReadBase(Base.C)
      .setReferenceBase(Base.C)
      .setMapQuality(30)
      .setSangerQuality(30)
      .setCountAtPosition(1)
      .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.C)
        .setReferenceBase(Base.C)
        .setMapQuality(40)
        .setSangerQuality(40)
        .setCountAtPosition(1)
        .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.A)
        .setReferenceBase(Base.C)
        .setMapQuality(40)
        .setSangerQuality(30)
        .setCountAtPosition(1)
        .build())
    //TODO(arahuja) test is not correct as multiplying probabilities is not exactly the same as averaging qual scores
    val expected = List((8.0 * ((0.999 * 0.999) * (0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999))) / 8.0,

      ((0.999 * 0.999 + (1.0 - 0.999 * 0.999)) * (0.9999 * 0.9999 + (1.0 - 0.9999 * 0.9999))
        * ((1.0 - 0.999 * 0.9999) + (0.999 * 0.9999))) / 8.0,
      (8.0 * (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (0.999 * 0.9999)) / 8.0).reverse

    val scored = call.scoreGenotypeLikelihoods(pl)

    for (i <- 0 to 2) {
      assertAlmostEqual(expected(i), scored(i), 1e-3)
    }
  }

  test("score genotype for single sample, all bases non-ref") {
    val call = new PileupCallSimpleSNP(2)

    val pl = List(ADAMPileup.newBuilder()
      .setReadBase(Base.A)
      .setReferenceBase(Base.C)
      .setMapQuality(30)
      .setSangerQuality(30)
      .setCountAtPosition(1)
      .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.A)
        .setReferenceBase(Base.C)
        .setMapQuality(40)
        .setSangerQuality(40)
        .setCountAtPosition(1)
        .build(),
      ADAMPileup.newBuilder()
        .setReadBase(Base.A)
        .setReferenceBase(Base.C)
        .setMapQuality(40)
        .setSangerQuality(30)
        .setCountAtPosition(1)
        .build())

    val expected = List(8.0 * (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999) / 8.0,
      ((0.999 * 0.999 * 0.9999 * 0.9999 * 0.999 * 0.9999) +
        (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999)) / 8.0,
      (8.0 * (0.999 * 0.999 * 0.9999 * 0.9999 * 0.999 * 0.9999)) / 8.0).reverse

    val scored = call.scoreGenotypeLikelihoods(pl)

    for (i <- 0 to 2) {
      assertAlmostEqual(expected(i), scored(i))
    }
  }

}
