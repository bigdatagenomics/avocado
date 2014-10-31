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
package org.bdgenomics.avocado.algorithms.reference

import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.formats.avro.{ Contig, NucleotideContigFragment }
import org.scalatest.FunSuite

class ResizeAndFlankReferenceFragmentsSuite extends FunSuite {

  test("resize a small fragment") {
    val fragments = ResizeAndFlankReferenceFragments.resize((ReferenceRegion("chr1", 0L, 10L),
      NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setFragmentSequence("AAAAATTTTT")
      .setFragmentStartPosition(0L)
      .build()), 5).map(_._2)

    assert(fragments.length === 2)
    fragments.foreach(f => assert(f.getFragmentSequence.length === 5))
    assert(ReferenceRegion(fragments(0)).get === ReferenceRegion("chr1", 0L, 5L))
    assert(fragments(0).getFragmentSequence === "AAAAA")
    assert(ReferenceRegion(fragments(1)).get === ReferenceRegion("chr1", 5L, 10L))
    assert(fragments(1).getFragmentSequence === "TTTTT")
  }

  test("don't put flanks on non-adjacent fragments") {
    val testIter = Iterator((ReferenceRegion("chr1", 0L, 10L),
      NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setFragmentSequence("AAAAATTTTT")
      .setFragmentStartPosition(0L)
      .build()), (ReferenceRegion("chr1", 20L, 30L),
      NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setFragmentSequence("CCCCCGGGGG")
      .setFragmentStartPosition(20L)
      .build()))

    val fragments = ResizeAndFlankReferenceFragments.overlapAndResize(testIter, 10, 5).toSeq

    assert(fragments.size === 2)
    fragments.foreach(_.getFragmentSequence.length === 10)
    assert(fragments(0).getFragmentSequence === "AAAAATTTTT")
    assert(fragments(0).getFragmentStartPosition === 0L)
    assert(fragments(1).getFragmentSequence === "CCCCCGGGGG")
    assert(fragments(1).getFragmentStartPosition === 20L)
  }

  test("put flanks on adjacent fragments") {
    val testIter = Iterator((ReferenceRegion("chr1", 0L, 10L),
      NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setFragmentSequence("AAAAATTTTT")
      .setFragmentStartPosition(0L)
      .build()), (ReferenceRegion("chr1", 10L, 20L),
      NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setFragmentSequence("NNNNNUUUUU")
      .setFragmentStartPosition(10L)
      .build()), (ReferenceRegion("chr1", 20L, 30L),
      NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setFragmentSequence("CCCCCGGGGG")
      .setFragmentStartPosition(20L)
      .build()))

    val fragments = ResizeAndFlankReferenceFragments.overlapAndResize(testIter, 10, 5).toSeq

    assert(fragments.size === 3)
    assert(fragments(0).getFragmentSequence === "AAAAATTTTTNNNNN")
    assert(fragments(0).getFragmentStartPosition === 0L)
    assert(fragments(1).getFragmentSequence === "TTTTTNNNNNUUUUUCCCCC")
    assert(fragments(1).getFragmentStartPosition === 5L)
    assert(fragments(2).getFragmentSequence === "UUUUUCCCCCGGGGG")
    assert(fragments(2).getFragmentStartPosition === 15L)
  }

  test("break and flank a single fragment") {
    val testIter = Iterator((ReferenceRegion("chr1", 0L, 10L),
      NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setFragmentSequence("AAAAATTTTT")
      .setFragmentStartPosition(0L)
      .build()))

    val fragments = ResizeAndFlankReferenceFragments.overlapAndResize(testIter, 5, 2).toSeq

    assert(fragments.size === 2)
    fragments.foreach(_.getFragmentSequence.length === 5)
    assert(fragments(0).getFragmentSequence === "AAAAATT")
    assert(fragments(0).getFragmentStartPosition === 0L)
    assert(fragments(1).getFragmentSequence === "AATTTTT")
    assert(fragments(1).getFragmentStartPosition === 3L)
  }
}
