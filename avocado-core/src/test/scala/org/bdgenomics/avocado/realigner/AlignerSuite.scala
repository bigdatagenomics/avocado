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
package org.bdgenomics.avocado.realigner

import org.bdgenomics.avocado.models.{
  Deletion,
  Insertion,
  Match
}
import org.scalatest.FunSuite

class AlignerSuite extends FunSuite {

  test("align a minimally flanked sequence with a snp") {

    // k = 10
    //                 |
    // read: ACTCTGAAACTAACCACTGGA
    // ref:  ACTCTGAAACAAACCACTGGA
    val aln = Aligner.align("ACTCTGAAACAAACCACTGGA",
      "ACTCTGAAACTAACCACTGGA",
      10)

    assert(aln.size === 3)
    assert(aln(0) === Match(10))
    assert(aln(1) === Match(1, Some("A")))
    assert(aln(2) === Match(10))
  }

  test("align a minimally flanked sequence with a 3 bp mnp") {

    // k = 10
    //                 |||
    // read: ACTCTGAAACTCCAACCACTGGA
    // ref:  ACTCTGAAACAGAAACCACTGGA
    val aln = Aligner.align("ACTCTGAAACAGAAACCACTGGA",
      "ACTCTGAAACTCCAACCACTGGA",
      10)

    assert(aln.size === 3)
    assert(aln(0) === Match(10))
    assert(aln(1) === Match(3, Some("AGA")))
    assert(aln(2) === Match(10))
  }

  test("align a minimally flanked sequence with 2 snps separated by 1bp") {

    // k = 10
    //                 | |
    // read: ACTCTGAAACTGCAACCACTGGA
    // ref:  ACTCTGAAACAGAAACCACTGGA
    val aln = Aligner.align("ACTCTGAAACAGAAACCACTGGA",
      "ACTCTGAAACTGCAACCACTGGA",
      10)

    assert(aln.size === 5)
    assert(aln(0) === Match(10))
    assert(aln(1) === Match(1, Some("A")))
    assert(aln(2) === Match(1))
    assert(aln(3) === Match(1, Some("A")))
    assert(aln(4) === Match(10))
  }

  test("align a minimally flanked sequence with 2 snps separated by 3bp") {

    // k = 10
    //                 |   |
    // read: ACTCTGAAACTGCCCAACCACTGGA
    // ref:  ACTCTGAAACAGCCAAACCACTGGA
    val aln = Aligner.align("ACTCTGAAACAGCCAAACCACTGGA",
      "ACTCTGAAACTGCCCAACCACTGGA",
      10)

    assert(aln.size === 5)
    assert(aln(0) === Match(10))
    assert(aln(1) === Match(1, Some("A")))
    assert(aln(2) === Match(3))
    assert(aln(3) === Match(1, Some("A")))
    assert(aln(4) === Match(10))
  }

  test("align a minimally flanked sequence with a simple insert") {

    // k = 10
    //                 
    // read: CTCTGAAACAACAACCACTGG
    // ref:  CTCTGAAACA__ACCACTGG
    val aln = Aligner.align("CTCTGAAACAAACCACTGG",
      "CTCTGAAACAACAACCACTGG",
      10)

    assert(aln.size === 3)
    assert(aln(0) === Match(10))
    assert(aln(1) === Insertion(2))
    assert(aln(2) === Match(9))
  }

  test("align a minimally flanked sequence with a complex insert") {

    // k = 10
    //                 
    // read: CTCTGAAACAACAACCACTGGT
    // ref:  CTCTGAAACA_TAACCACTGGT
    val aln = Aligner.align("CTCTGAAACATAACCACTGGT",
      "CTCTGAAACAACAACCACTGGT",
      10)

    assert(aln.size === 4)
    assert(aln(0) === Match(10))
    assert(aln(1) === Insertion(2))
    assert(aln(2) === Deletion("T"))
    assert(aln(3) === Match(10))
  }

  test("align a minimally flanked sequence with a simple deletion") {

    // k = 10
    //                 
    // read: ACTCTGAATA___AACCACTGGA
    // ref:  ACTCTGAATAGGGAACCACTGGA
    val aln = Aligner.align("ACTCTGAATAGGGAACCACTGGA",
      "ACTCTGAATAAACCACTGGA",
      10)

    assert(aln.size === 3)
    assert(aln(0) === Match(10))
    assert(aln(1) === Deletion("GGG"))
    assert(aln(2) === Match(10))
  }

  test("align a minimally flanked sequence that contains a discordant k-mer pair") {
    //         RPT  RPT
    // alt: ATGGGACAACCGAA
    // ref: ATGACCTGGGAGAA
    //         ||||||||
    val aln = Aligner.align("ATGGGACAACCGAA",
      "ATGACCTGGGAGAA",
      3)

    assert(aln.size === 3)
    assert(aln(0) === Match(3))
    assert(aln(1) === Match(8, Some("GGACAACC")))
    assert(aln(2) === Match(3))
  }

  test("align a minimally flanked sequence with a complex deletion") {

    // k = 10
    //                 
    // read: ACTCTGAATA__CAACCACTGGA
    // ref:  ACTCTGAATAGGGAACCACTGGA
    //                   |
    val aln = Aligner.align("ACTCTGAATAGGGAACCACTGGA",
      "ACTCTGAATACAACCACTGGA",
      10)

    assert(aln.size === 4)
    assert(aln(0) === Match(10))
    assert(aln(1) === Deletion("GGG"))
    assert(aln(2) === Insertion(1))
    assert(aln(3) === Match(10))
  }

  test("align a minimally flanked sequence with 2 snps separated by two matching k-mers") {

    // k = 10
    //                 |           |
    // read: ACTCTGAAACTGCCTAGCAACACAACCACTGGA
    // ref:  ACTCTGAAACAGCCTAGCAACAAAACCACTGGA
    val aln = Aligner.align("ACTCTGAAACAGCCTAGCAACAAAACCACTGGA",
      "ACTCTGAAACTGCCTAGCAACACAACCACTGGA",
      10)

    assert(aln.size === 5)
    assert(aln(0) === Match(10))
    assert(aln(1) === Match(1, Some("A")))
    assert(aln(2) === Match(11))
    assert(aln(3) === Match(1, Some("A")))
    assert(aln(4) === Match(10))
  }

  test("align a minimally flanked sequence with a snp and an indel separated by one matching k-mer") {

    // k = 10
    //                 |
    // read: ACTCTGAAACTGCCTAGCACC_AAACCACTGGA
    // ref:  ACTCTGAAACAGCCTAGCACCAAAACCACTGGA
    val aln = Aligner.align("ACTCTGAAACAGCCTAGCACCAAAACCACTGGA",
      "ACTCTGAAACTGCCTAGCACCAAACCACTGGA",
      10)

    assert(aln.size === 5)
    assert(aln(0) === Match(10))
    assert(aln(1) === Match(1, Some("A")))
    assert(aln(2) === Match(10))
    assert(aln(3) === Deletion("A"))
    assert(aln(4) === Match(11))
  }

  test("zip and trim short insert") {
    // ref: CTCTGAAACA__AACCACTGG
    // alt: CTCTGAAACAACAACCACTGG
    //      sssssssssssIIeeeeeeee
    val (ref, alt, trimStart, trimEnd) = Aligner.zipAndTrim(
      "CTCTGAAACAAACCACTGG", "CTCTGAAACAACAACCACTGG")

    assert(trimStart === 10)
    assert(trimEnd === 9)
    assert(ref.isEmpty)
    assert(alt === "AC")
  }

  test("zip and trim short deletion") {
    // ref: ACTCTGAATAGGGAACCACTGGA
    // alt: ACTCTGAATA___AACCACTGGA
    //      ssssssssssDDDeeeeeeeeee
    val (ref, alt, trimStart, trimEnd) = Aligner.zipAndTrim(
      "ACTCTGAATAGGGAACCACTGGA", "ACTCTGAATAAACCACTGGA")

    assert(trimStart === 10)
    assert(trimEnd === 10)
    assert(ref === "GGG")
    assert(alt.isEmpty)
  }

  test("cut up a sequence that is longer than the k-mer length") {
    val kmerMap = Aligner.toKmers("ACACTGTG", 4)

    assert(kmerMap.size === 5)
    assert(kmerMap("ACAC") === 0)
    assert(kmerMap("CACT") === 1)
    assert(kmerMap("ACTG") === 2)
    assert(kmerMap("CTGT") === 3)
    assert(kmerMap("TGTG") === 4)
  }

  test("cutting up a sequence that is shorter than the k-mer length yields an empty map") {
    val kmerMap = Aligner.toKmers("ACA", 4)

    assert(kmerMap.isEmpty)
  }

  test("cutting up a repeated sequence throws an assert") {
    intercept[AssertionError] {
      Aligner.toKmers("ACACTTACAC", 4)
    }
  }

  test("get no indices if we have no intersection") {
    // ref: ACCTG
    //      ACC
    //       CCT
    //        CTG
    // alt: GACTA
    //      GAC
    //       ACT
    //        CTA
    val refKmers = Aligner.toKmers("ACCTG", 3)
    val altKmers = Aligner.toKmers("GACTA", 3)

    val indices = Aligner.kmerIndices(refKmers.keySet & altKmers.keySet,
      refKmers,
      altKmers)

    assert(indices.isEmpty)
  }

  test("get correct index for a single intersection") {
    // ref: ACCTG
    //      ACC
    //       CCT   - match 1
    //        CTG
    // alt: GCCCTA
    //      GCC
    //       CCC
    //        CCT  - match 1
    //         CTA
    val refKmers = Aligner.toKmers("ACCTG", 3)
    val altKmers = Aligner.toKmers("GCCCTA", 3)

    val indices = Aligner.kmerIndices(refKmers.keySet & altKmers.keySet,
      refKmers,
      altKmers)

    assert(indices.size === 1)
    assert(indices(0) === (1, 2))
  }

  test("get correct indices for two k-mers in a row") {
    // ref: ACCTAG
    //      ACC
    //       CCT   - match 1
    //        CTA  - match 2
    //         TAG
    // alt: GCCCTA
    //      GCC
    //       CCC
    //        CCT  - match 1
    //         CTA - match 2
    val refKmers = Aligner.toKmers("ACCTAG", 3)
    val altKmers = Aligner.toKmers("GCCCTA", 3)

    val indices = Aligner.kmerIndices(refKmers.keySet & altKmers.keySet,
      refKmers,
      altKmers)

    assert(indices.size === 2)
    assert(indices(0) === (1, 2))
    assert(indices(1) === (2, 3))
  }

  test("get correct indices for two k-mers separated by a snp") {
    // ref: ACCTGGCA
    //      ACC
    //       CCT      - match 1
    //        CTG
    //         TGG
    //          GGC
    //           GCA  - match 2
    // alt: GCCCTAGCA
    //      GCC
    //       CCC
    //        CCT     - match 1
    //         CTA
    //          TAG
    //           AGC
    //            GCA - match 2
    val refKmers = Aligner.toKmers("ACCTGGCA", 3)
    val altKmers = Aligner.toKmers("GCCCTAGCA", 3)

    val indices = Aligner.kmerIndices(refKmers.keySet & altKmers.keySet,
      refKmers,
      altKmers)

    assert(indices.size === 2)
    assert(indices(0) === (1, 2))
    assert(indices(1) === (5, 6))
  }

  test("get correct indices for two k-mers separated by an indel") {
    // ref: ACCTGGCA
    //      ACC
    //       CCT       - match 1
    //        CTG      - match 2
    //         TGG
    //          GGC
    //           GCA   - match 3
    // alt: GCCCTGAGCA
    //      GCC
    //       CCC
    //        CCT      - match 1
    //         CTG     - match 2
    //          TGA
    //           GAG
    //            AGC
    //             GCA - match 3
    val refKmers = Aligner.toKmers("ACCTGGCA", 3)
    val altKmers = Aligner.toKmers("GCCCTGAGCA", 3)

    val indices = Aligner.kmerIndices(refKmers.keySet & altKmers.keySet,
      refKmers,
      altKmers)

    assert(indices.size === 3)
    assert(indices(0) === (1, 2))
    assert(indices(1) === (2, 3))
    assert(indices(2) === (5, 7))
  }

  test("get correct indices for two k-mers whose positions are flipped") {
    // ref: ACCTGGGC
    //      ACC      - match 1
    //       CCT
    //        CTG
    //         TGG
    //          GGG
    //           GGC - match 2
    // alt: GGCCAACC
    //      GGC      - match 2
    //       GCC
    //        CCA
    //         CAA
    //          AAC
    //           ACC - match 1
    val refKmers = Aligner.toKmers("ACCTGGGC", 3)
    val altKmers = Aligner.toKmers("GGCCAACC", 3)

    val indices = Aligner.kmerIndices(refKmers.keySet & altKmers.keySet,
      refKmers,
      altKmers)

    assert(indices.size === 2)
    assert(indices(0) === (0, 5))
    assert(indices(1) === (5, 0))
  }

  test("fire assert when checking negative index pair") {
    intercept[AssertionError] {
      val indexSeq = Seq((-1, -1))

      Aligner.indicesHaveConcordantOrder(indexSeq)
    }
  }

  test("a set of a single index pair is concordant") {
    val indexSeq = Seq((3, 6))

    assert(Aligner.indicesHaveConcordantOrder(indexSeq))
  }

  test("a set with a pair of index pairs is concordant") {
    val indexSeq = Seq((1, 2), (2, 3))

    assert(Aligner.indicesHaveConcordantOrder(indexSeq))
  }

  test("a set with multiple good index pairs is concordant") {
    val indexSeq = Seq((1, 2), (2, 3), (5, 7))

    assert(Aligner.indicesHaveConcordantOrder(indexSeq))
  }

  test("a set with a pair of swapped index pairs is discordant") {
    val indexSeq = Seq((0, 5), (5, 0))

    assert(!Aligner.indicesHaveConcordantOrder(indexSeq))
  }

  test("a set with a pair of both con/discordant index pairs is discordant") {
    val indexSeq = Seq((0, 8), (1, 9), (5, 0))

    assert(!Aligner.indicesHaveConcordantOrder(indexSeq))
  }

  test("making blocks from no indices returns a single unknown block") {
    val blocks = Aligner.indicesToBlocks(Seq.empty,
      "ACAT",
      "TCC",
      4)

    assert(blocks.size === 1)
    assert(blocks(0) === UnknownBlock("ACAT", "TCC"))
  }

  test("make blocks from a single match between two snps") {

    // ref: ACATAC
    // alt: GCATAG
    //      |    |
    val blocks = Aligner.indicesToBlocks(Seq((1, 1)),
      "ACATAC",
      "GCATAG",
      4)

    assert(blocks.size === 3)
    assert(blocks(0) === UnknownBlock("A", "G"))
    assert(blocks(1) === MatchBlock(4))
    assert(blocks(2) === UnknownBlock("C", "G"))
  }

  test("make blocks from three matches between two snps") {
    // ref: ACATACCC
    // alt: GCATACCG
    //      |      |
    val blocks = Aligner.indicesToBlocks(Seq((1, 1),
      (2, 2),
      (3, 3)),
      "ACATACCC",
      "GCATACCG",
      4)

    assert(blocks.size === 3)
    assert(blocks(0) === UnknownBlock("A", "G"))
    assert(blocks(1) === MatchBlock(6))
    assert(blocks(2) === UnknownBlock("C", "G"))
  }

  test("make blocks from three matches between two indels, opposite events") {
    // ref: ACATACC_
    // alt: _CATACCGT
    val blocks = Aligner.indicesToBlocks(Seq((1, 0),
      (2, 1),
      (3, 2)),
      "ACATACC",
      "CATACCGT",
      4)

    assert(blocks.size === 3)
    assert(blocks(0) === UnknownBlock("A", ""))
    assert(blocks(1) === MatchBlock(6))
    assert(blocks(2) === UnknownBlock("", "GT"))
  }

  test("make blocks from three matches between two indels, same events") {
    // ref: _CATAC__
    // alt: GCATACGT
    val blocks = Aligner.indicesToBlocks(Seq((0, 1),
      (1, 2)),
      "CATAC",
      "GCATACGT",
      4)

    assert(blocks.size === 3)
    assert(blocks(0) === UnknownBlock("", "G"))
    assert(blocks(1) === MatchBlock(5))
    assert(blocks(2) === UnknownBlock("", "GT"))
  }

  test("make blocks from matches between snp/indel/snp") {

    // ref: ACAT__ACTCGC <-- note: CATA k-mer matches pre-justification
    // alt: GCATAAACTCGG
    //      |          |
    val blocks = Aligner.indicesToBlocks(Seq((1, 1), (4, 6), (5, 7)),
      "ACATACTCGC",
      "GCATAAACTCGG",
      4)

    assert(blocks.size === 5)
    assert(blocks(0) === UnknownBlock("A", "G"))
    assert(blocks(1) === MatchBlock(3))
    assert(blocks(2) === UnknownBlock("", "AA"))
    assert(blocks(3) === MatchBlock(5))
    assert(blocks(4) === UnknownBlock("C", "G"))
  }
}
