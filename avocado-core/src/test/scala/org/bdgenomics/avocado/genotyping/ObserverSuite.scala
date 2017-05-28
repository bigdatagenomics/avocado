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

import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.utils.misc.MathUtils
import scala.math.{ log, pow, sqrt }

class ObserverSuite extends AvocadoFunSuite {

  lazy val summaryObservations = ScoredObservation.createScores(sc, 90, 90, 2)
    .collect()
    .toSeq

  test("a fully clipped read will not generate any observations") {
    val read = AlignmentRecord.newBuilder
      .setStart(10L)
      .setEnd(11L)
      .setContigName("ctg")
      .setSequence("AAAA")
      .setQual("****")
      .setMapq(0)
      .setReadMapped(true)
      .setCigar("4S")
      .setMismatchingPositions("0")
      .build()

    assert(Observer.observeRead(read, 2).isEmpty)
  }

  def eps(g: Int,
          mapP: Double,
          baseP: Double): Double = {
    g * (1.0 - (mapP * baseP))
  }

  def invEps(g: Int,
             mapP: Double,
             baseP: Double): Double = {
    g * mapP * baseP
  }

  def logL(g: Int,
           m: Int,
           mapP: Double,
           baseP: Double): Double = {
    log(eps(m - g, mapP, baseP) + invEps(g, mapP, baseP))
  }

  sparkTest("generate observations for a sequence match under diploid model") {
    val read = AlignmentRecord.newBuilder
      .setStart(10L)
      .setEnd(15L)
      .setContigName("ctg")
      .setSequence("ACGT")
      .setQual(Array(20, 30, 40, 50)
        .map(v => (v + 33).toChar)
        .mkString)
      .setReadMapped(true)
      .setReadNegativeStrand(false)
      .setCigar("4M")
      .setMismatchingPositions("4")
      .setMapq(50)
      .setRecordGroupSample("sample")
      .build()

    val obs = Observer.observeRead(read, 2)
      .toSeq
      .map(p => (p._1, p._2.toObservation(summaryObservations)))

    assert(obs.size === 4)
    assert(obs.forall(_._2.alleleCoverage == 0))
    assert(obs.forall(_._2.otherCoverage == 1))
    assert(obs.forall(_._2.alleleForwardStrand == 0))
    assert(obs.forall(_._2.otherForwardStrand == 1))
    assert(obs.forall(_._2.alleleLogLikelihoods.size == 3))
    assert(obs.forall(_._2.otherLogLikelihoods.size == 3))
    assert(obs.forall(o => {
      MathUtils.fpEquals(o._2.squareMapQ, pow(50.0, 2.0))
    }))
    assert(obs.forall(_._2.isRef))
    assert(obs.forall(_._1._3 == "sample"))
    (0 until 4).zip("ACGT").foreach(p => {
      val (idx, base) = p
      assert(obs(idx)._1._1 === ReferenceRegion("ctg", 10 + idx, 11 + idx))
      assert(obs(idx)._1._2 === base.toString)
    })
    (0 until 4).zip(Array(0.99, 0.999, 0.9999, 0.99999)).foreach(p => {
      val (idx, baseQ) = p
      assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
    })
  }

  test("generate observations for a read with an insertion under diploid model") {
    val read = AlignmentRecord.newBuilder
      .setStart(10L)
      .setEnd(12L)
      .setContigName("ctg")
      .setSequence("ACGT")
      .setQual(Array(20, 30, 40, 50)
        .map(v => (v + 33).toChar)
        .mkString)
      .setReadMapped(true)
      .setReadNegativeStrand(false)
      .setCigar("1M2I1M")
      .setMismatchingPositions("2")
      .setMapq(50)
      .setRecordGroupSample("sample")
      .build()

    val obs = Observer.observeRead(read, 2)
      .toSeq
      .map(p => (p._1, p._2.toObservation(summaryObservations)))

    assert(obs.size === 3)
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.alleleCoverage == 1))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.alleleCoverage == 0))
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.otherCoverage == 0))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.otherCoverage == 1))
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.alleleForwardStrand == 1))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.alleleForwardStrand == 0))
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.otherForwardStrand == 0))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.otherForwardStrand == 1))
    assert(obs.forall(_._2.alleleLogLikelihoods.size == 3))
    assert(obs.forall(_._2.otherLogLikelihoods.size == 3))
    assert(obs.count(_._2.isRef) === 2)
    assert(obs.forall(o => {
      MathUtils.fpEquals(o._2.squareMapQ, pow(50.0, 2.0))
    }))
    assert(obs.forall(_._1._3 == "sample"))
    Seq((0, 0), (1, 2)).zip("AT").foreach(p => {
      val ((offset, idx), base) = p
      assert(obs(idx)._1._1 === ReferenceRegion("ctg", 10 + offset, 11 + offset))
      assert(obs(idx)._1._2 === base.toString)
    })
    assert(obs(1)._1._1 === ReferenceRegion("ctg", 10, 11))
    assert(obs(1)._1._2 === "CG")
    (0 until 3).zip(Array(0.99, 1.0 - sqrt(0.001 * 0.0001), 0.99999)).foreach(p => {
      val (idx, baseQ) = p
      if (obs(idx)._1._2.size != 1) {
        assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
      } else {
        assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
      }
    })
  }

  test("generate observations for a read with a deletion under diploid model") {
    val read = AlignmentRecord.newBuilder
      .setStart(10L)
      .setEnd(17L)
      .setContigName("ctg")
      .setSequence("ACGT")
      .setQual(Array(20, 30, 40, 50)
        .map(v => (v + 33).toChar)
        .mkString)
      .setReadMapped(true)
      .setReadNegativeStrand(false)
      .setCigar("2M2D2M")
      .setMismatchingPositions("2^NN2")
      .setMapq(50)
      .setRecordGroupSample("sample")
      .build()

    val obs = Observer.observeRead(read, 2)
      .toSeq
      .map(p => (p._1, p._2.toObservation(summaryObservations)))

    assert(obs.size === 5)
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.alleleCoverage == 1))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.alleleCoverage == 0))
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.otherCoverage == 0))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.otherCoverage == 1))
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.alleleForwardStrand == 1))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.alleleForwardStrand == 0))
    assert(obs.filter(p => p._1._2.length != 1).forall(_._2.otherForwardStrand == 0))
    assert(obs.filter(p => p._1._2.length == 1).forall(_._2.otherForwardStrand == 1))
    assert(obs.forall(_._2.alleleLogLikelihoods.size == 3))
    assert(obs.forall(_._2.otherLogLikelihoods.size == 3))
    assert(obs.count(_._2.isRef) === 4)
    assert(obs.forall(o => {
      MathUtils.fpEquals(o._2.squareMapQ, pow(50.0, 2.0))
    }))
    assert(obs.forall(_._1._3 == "sample"))
    Seq((0, 0), (1, 1), (3, 4), (4, 5)).zip("ACGT").foreach(p => {
      val ((idx, pos), base) = p
      assert(obs(idx)._1._1 === ReferenceRegion("ctg", 10 + pos, 11 + pos))
      assert(obs(idx)._1._2 === base.toString)
    })
    assert(obs(2)._1._1 === ReferenceRegion("ctg", 12, 14))
    assert(obs(2)._1._2.isEmpty)
    (0 until 5).zip(Array(0.99, 0.999, 1.0, 0.9999, 0.99999)).foreach(p => {
      val (idx, baseQ) = p
      if (obs(idx)._1._2.size != 1) {
        assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
      } else {
        assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
        assert(MathUtils.fpEquals(obs(idx)._2.referenceLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
      }
    })
  }
}
