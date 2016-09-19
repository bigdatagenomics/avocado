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
package org.bdgenomics.avocado.models

import org.scalatest.FunSuite

class ObservationOperatorSuite extends FunSuite {

  test("zero operators are empty") {
    assert(!Match(0).nonEmpty)
    assert(!Insertion(0).nonEmpty)
    assert(!Deletion("").nonEmpty)
    assert(!Clipped(0).nonEmpty)
  }

  test("non-zero operators are non-empty") {
    assert(Match(2).nonEmpty)
    assert(Insertion(2).nonEmpty)
    assert(Deletion("AB").nonEmpty)
    assert(Clipped(4).nonEmpty)
  }

  test("cannot build mismatch with wrong ref length") {
    intercept[IllegalArgumentException] {
      Match(1, Some("AC"))
    }
  }

  test("collapsing a non repeated set of operators should eliminate 0 ops") {
    val obs = ObservationOperator.collapse(Iterable(
      Match(10),
      Insertion(5),
      Deletion(""),
      Match(10),
      Insertion(3),
      Match(0),
      Deletion("AACAG"),
      Match(10))).toSeq

    assert(obs.size === 6)
    assert(obs(0) === Match(10))
    assert(obs(1) === Insertion(5))
    assert(obs(2) === Match(10))
    assert(obs(3) === Insertion(3))
    assert(obs(4) === Deletion("AACAG"))
    assert(obs(5) === Match(10))
    assert(obs.mkString === "10=5I10=3I5D10=")
  }

  test("collapsing a repeated set of operators with mixed match/mismatch") {
    val obs = ObservationOperator.collapse(Iterable(
      Match(10),
      Insertion(5),
      Match(5),
      Match(1, Some("A")),
      Match(1, Some("C")),
      Match(3),
      Insertion(3),
      Match(0),
      Deletion("AACAG"),
      Match(10))).toSeq

    assert(obs.size === 8)
    assert(obs(0) === Match(10))
    assert(obs(1) === Insertion(5))
    assert(obs(2) === Match(5))
    assert(obs(3) === Match(2, Some("AC")))
    assert(obs(4) === Match(3))
    assert(obs(5) === Insertion(3))
    assert(obs(6) === Deletion("AACAG"))
    assert(obs(7) === Match(10))
    assert(obs.mkString === "10=5I5=2X3=3I5D10=")
  }

  test("collapse a set of operators with repeats") {
    val obs = ObservationOperator.collapse(Iterable(
      Match(10),
      Insertion(5),
      Deletion(""),
      Insertion(5),
      Match(10),
      Match(5),
      Insertion(3),
      Deletion("GAG"),
      Match(0),
      Deletion(""),
      Match(0),
      Deletion("ACACG"),
      Deletion("TC"),
      Match(10))).toSeq

    assert(obs.size === 6)
    assert(obs(0) === Match(10))
    assert(obs(1) === Insertion(10))
    assert(obs(2) === Match(15))
    assert(obs(3) === Insertion(3))
    assert(obs(4) === Deletion("GAGACACGTC"))
    assert(obs(5) === Match(10))
    assert(obs.mkString === "10=10I15=3I10D10=")
  }

  test("collapse a set of operators with repeats and clips") {
    val obs = ObservationOperator.collapse(Iterable(
      Clipped(3, soft = false),
      Clipped(5),
      Match(10),
      Insertion(5),
      Deletion(""),
      Insertion(5),
      Match(10),
      Match(5),
      Insertion(3),
      Deletion("GAG"),
      Match(0),
      Deletion(""),
      Match(0),
      Deletion("ACACG"),
      Deletion("TC"),
      Match(10),
      Clipped(3))).toSeq

    assert(obs.size === 9)
    assert(obs(0) === Clipped(3, soft = false))
    assert(obs(1) === Clipped(5))
    assert(obs(2) === Match(10))
    assert(obs(3) === Insertion(10))
    assert(obs(4) === Match(15))
    assert(obs(5) === Insertion(3))
    assert(obs(6) === Deletion("GAGACACGTC"))
    assert(obs(7) === Match(10))
    assert(obs(8) === Clipped(3))
    assert(obs.mkString === "3H5S10=10I15=3I10D10=3S")
  }

  test("make a cigar and md tag from a single sequence match") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Match(10)))

    assert(cigar === "10=")
    assert(md === "10")
  }

  test("make a cigar and md tag from a single sequence mismatch") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Match(1, Some("A"))))

    assert(cigar === "1X")
    assert(md === "0A0")
  }

  test("make a cigar and md tag from a single multi-base sequence match") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Match(2, Some("AT"))))

    assert(cigar === "2X")
    assert(md === "0A0T0")
  }

  test("make a cigar and md tag from a single deletion") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Deletion("A")))

    assert(cigar === "1D")
    assert(md === "0^A0")
  }

  test("make a cigar and md tag from a single insertion") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Insertion(4)))

    assert(cigar === "4I")
    assert(md === "0")
  }

  test("make a cigar for a match followed by a deletion") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Match(4),
      Deletion("AC")))

    assert(cigar === "4=2D")
    assert(md === "4^AC0")
  }

  test("make a cigar for an insertion flanked by matches") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Match(4),
      Insertion(4),
      Match(6)))

    assert(cigar === "4=4I6=")
    assert(md === "10")
  }

  test("make a cigar for a match followed by a mismatch") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Match(4),
      Match(1, Some("T"))))

    assert(cigar === "4=1X")
    assert(md === "4T0")
  }

  test("make a cigar for a multi-base mismatch flanked by matches") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Match(8),
      Match(2, Some("TC")),
      Match(2)))

    assert(cigar === "8=2X2=")
    assert(md === "8T0C2")
  }

  test("make a cigar for a match after a clip") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Clipped(2),
      Match(8)))

    assert(cigar === "2S8=")
    assert(md === "8")
  }

  test("make a cigar for a mismatch after a clip") {
    val (cigar, md) = ObservationOperator.makeCigarAndMD(Iterable(Clipped(1, soft = false),
      Match(1, Some("G"))))

    assert(cigar === "1H1X")
    assert(md === "0G0")
  }

  test("extract reference from a single snp") {
    val ref = ObservationOperator.extractReference("A",
      Iterable(Match(1, Some("T"))))

    assert(ref === "T")
  }

  test("extract reference from a single deletion") {
    val ref = ObservationOperator.extractReference("",
      Iterable(Deletion("AT")))

    assert(ref === "AT")
  }

  test("extract reference from a single insertion") {
    val ref = ObservationOperator.extractReference("ACAG",
      Iterable(Match(1),
        Insertion(2),
        Match(1)))

    assert(ref === "AG")
  }

  test("extract reference from a soft clipped sequence") {
    val ref = ObservationOperator.extractReference("ACAGT",
      Iterable(Clipped(2),
        Match(3)))

    assert(ref === "AGT")
  }

  test("extract reference from a hard clipped sequence") {
    val ref = ObservationOperator.extractReference("ACAGT",
      Iterable(Match(5),
        Clipped(2, soft = false)))

    assert(ref === "ACAGT")
  }

  test("extract reference from a match flanked deletion") {
    val ref = ObservationOperator.extractReference("TCAACAGT",
      Iterable(Match(3),
        Deletion("ACA"),
        Match(5),
        Clipped(2, soft = false)))
    assert(ref === "TCAACAACAGT")
  }

  test("extract reference from a match flanked insertion") {
    val ref = ObservationOperator.extractReference("TCAACAGCAGT",
      Iterable(Clipped(3),
        Match(3),
        Insertion(3),
        Match(2)))
    assert(ref === "ACAGT")
  }
}
