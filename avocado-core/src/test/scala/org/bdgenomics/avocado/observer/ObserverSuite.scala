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
package org.bdgenomics.avocado.observer

import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.utils.misc.MathUtils
import scala.math.{ log, pow, sqrt }

class ObserverSuite extends AvocadoFunSuite {

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

  test("generate observations for a sequence match under diploid model") {
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

    assert(obs.size === 4)
    assert(obs.forall(_._2.coverage == 1))
    assert(obs.forall(_._2.forwardStrand == 1))
    assert(obs.forall(_._2.alleleLogLikelihoods.size == 3))
    assert(obs.forall(_._2.otherLogLikelihoods.size == 3))
    assert(obs.forall(o => {
      MathUtils.fpEquals(o._2.squareMapQ, pow(50.0, 2.0))
    }))
    assert(obs.forall(_._1._3 == "sample"))
    (0 until 4).zip("ACGT").foreach(p => {
      val (idx, base) = p
      assert(obs(idx)._1._1 === ReferenceRegion("ctg", 10 + idx, 11 + idx))
      assert(obs(idx)._1._2 === base.toString)
    })
    (0 until 4).zip(Array(0.99, 0.999, 0.9999, 0.99999)).foreach(p => {
      val (idx, baseQ) = p
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
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

    assert(obs.size === 3)
    assert(obs.forall(_._2.coverage == 1))
    assert(obs.forall(_._2.forwardStrand == 1))
    assert(obs.forall(_._2.alleleLogLikelihoods.size == 3))
    assert(obs.forall(_._2.otherLogLikelihoods.size == 3))
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
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
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

    assert(obs.size === 5)
    assert(obs.forall(_._2.coverage == 1))
    assert(obs.forall(_._2.forwardStrand == 1))
    assert(obs.forall(_._2.alleleLogLikelihoods.size == 3))
    assert(obs.forall(_._2.otherLogLikelihoods.size == 3))
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
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(0), logL(0, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(1), logL(1, 2, baseQ, 0.99999)))
      assert(MathUtils.fpEquals(obs(idx)._2.alleleLogLikelihoods(2), logL(2, 2, baseQ, 0.99999)))
    })
  }

  sparkTest("unmapped reads generate no observations") {
    val reads = (0 until 100).map(i => {
      AlignmentRecord.newBuilder
        .setReadName(i.toString)
        .setReadMapped(false)
        .build()
    })

    val obsRdd = Observer.observe(sc.parallelize(reads), 2)
    assert(obsRdd.count === 0)
  }

  sparkTest("generate reads covering a sequence match") {
    val reads = Seq(AlignmentRecord.newBuilder
      .setRecordGroupSample("sample")
      .setStart(10L)
      .setEnd(18L)
      .setContigName("ctg")
      .setSequence("ACGTACGT")
      .setQual((20 + 33).toChar.toString * 8)
      .setReadMapped(true)
      .setReadNegativeStrand(false)
      .setCigar("8M")
      .setMismatchingPositions("8")
      .setMapq(40)
      .build(), AlignmentRecord.newBuilder
      .setRecordGroupSample("sample")
      .setStart(11L)
      .setEnd(19L)
      .setContigName("ctg")
      .setSequence("CGTACGTA")
      .setQual((30 + 33).toChar.toString * 8)
      .setReadMapped(true)
      .setReadNegativeStrand(true)
      .setCigar("8M")
      .setMismatchingPositions("8")
      .setMapq(50)
      .build(), AlignmentRecord.newBuilder
      .setRecordGroupSample("sample")
      .setStart(12L)
      .setEnd(20L)
      .setContigName("ctg")
      .setSequence("GTACGTAC")
      .setQual((40 + 33).toChar.toString * 8)
      .setReadMapped(true)
      .setReadNegativeStrand(false)
      .setCigar("8M")
      .setMismatchingPositions("8")
      .setMapq(60)
      .build())

    val r0l0 = logL(0, 2, 0.99, 0.9999)
    val r0l1 = logL(1, 2, 0.99, 0.9999)
    val r0l2 = logL(2, 2, 0.99, 0.9999)
    val r1l0 = logL(0, 2, 0.999, 0.99999)
    val r1l1 = logL(1, 2, 0.999, 0.99999)
    val r1l2 = logL(2, 2, 0.999, 0.99999)
    val r2l0 = logL(0, 2, 0.9999, 0.999999)
    val r2l1 = logL(1, 2, 0.9999, 0.999999)
    val r2l2 = logL(2, 2, 0.9999, 0.999999)

    val obsRdd = Observer.observe(sc.parallelize(reads), 2)
      .collect
    assert(obsRdd.length === 10)

    val tripleBases = obsRdd.filter(o => {
      val rr = o._1._1
      rr.start >= 12L && rr.end <= 18L
    })
    assert(tripleBases.length === 6)
    tripleBases.foreach(o => {
      val ((rr, allele, sample), obs) = o
      assert(rr.referenceName === "ctg")
      assert(sample === "sample")
      assert(allele.length === 1)
      assert(rr.length === 1)
      assert(obs.squareMapQ === 40 * 40 + 50 * 50 + 60 * 60)
      assert(obs.coverage === 3)
      assert(obs.forwardStrand === 2)
      assert(obs.alleleLogLikelihoods.length === 3)
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(0), r0l0 + r1l0 + r2l0))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(1), r0l1 + r1l1 + r2l1))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(2), r0l2 + r1l2 + r2l2))
    })

    val only1Bases = obsRdd.filter(o => {
      val rr = o._1._1
      rr.start == 10L && rr.end == 11L
    })
    assert(only1Bases.size === 1)
    only1Bases.foreach(o => {
      val ((rr, allele, sample), obs) = o
      assert(rr.referenceName === "ctg")
      assert(sample === "sample")
      assert(allele.length === 1)
      assert(rr.length === 1)
      assert(obs.squareMapQ === 40 * 40)
      assert(obs.coverage === 1)
      assert(obs.forwardStrand === 1)
      assert(obs.alleleLogLikelihoods.length === 3)
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(0), r0l0))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(1), r0l1))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(2), r0l2))
    })

    val both12Bases = obsRdd.filter(o => {
      val rr = o._1._1
      rr.start == 11L && rr.end == 12L
    })
    assert(both12Bases.size === 1)
    both12Bases.foreach(o => {
      val ((rr, allele, sample), obs) = o
      assert(rr.referenceName === "ctg")
      assert(sample === "sample")
      assert(allele.length === 1)
      assert(rr.length === 1)
      assert(obs.squareMapQ === 40 * 40 + 50 * 50)
      assert(obs.coverage === 2)
      assert(obs.forwardStrand === 1)
      assert(obs.alleleLogLikelihoods.length === 3)
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(0), r0l0 + r1l0))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(1), r0l1 + r1l1))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(2), r0l2 + r1l2))
    })

    val both23Bases = obsRdd.filter(o => {
      val rr = o._1._1
      rr.start == 18L && rr.end == 19L
    })
    assert(both23Bases.size === 1)
    both23Bases.foreach(o => {
      val ((rr, allele, sample), obs) = o
      assert(rr.referenceName === "ctg")
      assert(sample === "sample")
      assert(allele.length === 1)
      assert(rr.length === 1)
      assert(obs.squareMapQ === 50 * 50 + 60 * 60)
      assert(obs.coverage === 2)
      assert(obs.forwardStrand === 1)
      assert(obs.alleleLogLikelihoods.length === 3)
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(0), r2l0 + r1l0))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(1), r2l1 + r1l1))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(2), r2l2 + r1l2))
    })

    val only3Bases = obsRdd.filter(o => {
      val rr = o._1._1
      rr.start == 19L && rr.end == 20L
    })
    assert(only3Bases.size === 1)
    only3Bases.foreach(o => {
      val ((rr, allele, sample), obs) = o
      assert(rr.referenceName === "ctg")
      assert(sample === "sample")
      assert(allele.length === 1)
      assert(rr.length === 1)
      assert(obs.squareMapQ === 60 * 60)
      assert(obs.coverage === 1)
      assert(obs.forwardStrand === 1)
      assert(obs.alleleLogLikelihoods.length === 3)
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(0), r2l0))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(1), r2l1))
      assert(MathUtils.fpEquals(obs.alleleLogLikelihoods(2), r2l2))
    })
  }
}
