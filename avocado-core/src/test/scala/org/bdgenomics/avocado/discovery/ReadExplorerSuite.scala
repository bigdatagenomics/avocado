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
package org.bdgenomics.avocado.discovery

import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.formats.avro.{ AlignmentRecord, Contig }

class ReadExplorerSuite extends AvocadoFunSuite {

  sparkTest("observe a simple read") {
    val re = new ReadExplorer(sc.parallelize(Seq[Observation]()))

    val read = AlignmentRecord.newBuilder()
      .setStart(10L)
      .setEnd(15L)
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setMapq(40)
      .setSequence("ACTGA")
      .setQual("::;?:")
      .setCigar("5M")
      .setMismatchingPositions("2AA1")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations((read, 0L))
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })
    assert(observations.length === 5)
    assert(observations.forall(_.mapq == Some(40)))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.forall(_.distanceToNearestReadDeletion.isEmpty))
    assert(observations.forall(_.distanceToNearestReadInsertion.isEmpty))
    assert(observations.forall(_.alignedReadLen == 5))
    assert(observations.head.mismatchQScoreSum === Some(30 + 26))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 10L).head.phred === 25)
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.allele === "C")
    assert(observations.filter(_.pos.pos == 11L).head.phred === 25)
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.allele === "T")
    assert(observations.filter(_.pos.pos == 12L).head.phred === 26)
    assert(observations.filter(_.pos.pos == 13L).length === 1)
    assert(observations.filter(_.pos.pos == 13L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 13L).head.phred === 30)
    assert(observations.filter(_.pos.pos == 14L).length === 1)
    assert(observations.filter(_.pos.pos == 14L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 14L).head.phred === 25)

  }

  sparkTest("observe a read with a deletion") {
    val re = new ReadExplorer(sc.parallelize(Seq[Observation]()))

    val read = AlignmentRecord.newBuilder()
      .setStart(10L)
      .setEnd(17L)
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setMapq(40)
      .setSequence("ACTGA")
      .setQual(":::?:")
      .setCigar("2M2D3M")
      .setMismatchingPositions("2^AC1A1")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations((read, 0L))
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 5)
    assert(observations.filter(_.pos.pos != 15L).filter(_.allele != "_").forall(_.phred == 25))
    assert(observations.filter(_.allele == "_").forall(_.phred == 40))
    assert(observations.forall(_.mapq == Some(40)))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.exists(_.distanceToNearestReadDeletion.isDefined))
    assert(observations.exists(_.distanceToNearestReadInsertion.isEmpty))
    assert(observations.head.mismatchQScoreSum === Some(30))
    assert(observations.forall(_.alignedReadLen == 5))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.distanceToNearestReadDeletion === Some(2))
    assert(observations.filter(_.pos.pos == 10L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.allele === "C")
    assert(observations.filter(_.pos.pos == 11L).head.distanceToNearestReadDeletion === Some(1))
    assert(observations.filter(_.pos.pos == 11L).head.length === 3)
    assert(observations.filter(_.pos.pos == 12L).length === 0)
    assert(observations.filter(_.pos.pos == 13L).length === 0)
    assert(observations.filter(_.pos.pos == 14L).length === 1)
    assert(observations.filter(_.pos.pos == 14L).head.allele === "T")
    assert(observations.filter(_.pos.pos == 14L).head.distanceToNearestReadDeletion === Some(-1))
    assert(observations.filter(_.pos.pos == 15L).length === 1)
    assert(observations.filter(_.pos.pos == 15L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 15L).head.distanceToNearestReadDeletion === Some(-2))
    assert(observations.filter(_.pos.pos == 16L).length === 1)
    assert(observations.filter(_.pos.pos == 16L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 16L).head.distanceToNearestReadDeletion === Some(-3))

  }

  sparkTest("observe a read with an insertion") {
    val re = new ReadExplorer(sc.parallelize(Seq[Observation]()))

    val read = AlignmentRecord.newBuilder()
      .setStart(10L)
      .setEnd(13L)
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setMapq(40)
      .setSequence("ACTGA")
      .setQual(":::?:")
      .setCigar("1M2I2M")
      .setMismatchingPositions("1C1")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations((read, 0L))
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 3)
    assert(observations.filter(_.pos.pos != 11L).forall(_.phred == 25))
    assert(observations.forall(_.mapq == Some(40)))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.exists(_.distanceToNearestReadDeletion.isEmpty))
    assert(observations.exists(_.distanceToNearestReadInsertion.isDefined))
    assert(observations.forall(_.alignedReadLen == 5))
    assert(observations.head.mismatchQScoreSum === Some(30))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.distanceToNearestReadInsertion === Some(1))
    assert(observations.filter(_.pos.pos == 10L).head.allele === "ACT")
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.distanceToNearestReadInsertion === Some(-1))
    assert(observations.filter(_.pos.pos == 11L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 11L).head.phred === 30)
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.distanceToNearestReadInsertion === Some(-2))
    assert(observations.filter(_.pos.pos == 12L).head.allele === "A")
  }

  sparkTest("observe a simple read without map quality") {
    val re = new ReadExplorer(sc.parallelize(Seq[Observation]()))

    val read = AlignmentRecord.newBuilder()
      .setStart(10L)
      .setEnd(15L)
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setSequence("ACTGA")
      .setQual(":::::")
      .setCigar("5M")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations((read, 0L))
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 5)
    assert(observations.forall(_.phred == 25))
    assert(observations.forall(_.mapq == None))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.exists(_.distanceToNearestReadDeletion.isEmpty))
    assert(observations.exists(_.distanceToNearestReadInsertion.isEmpty))
    assert(observations.forall(_.alignedReadLen == 5))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.allele === "C")
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.allele === "T")
    assert(observations.filter(_.pos.pos == 13L).length === 1)
    assert(observations.filter(_.pos.pos == 13L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 14L).length === 1)
    assert(observations.filter(_.pos.pos == 14L).head.allele === "A")
  }

  sparkTest("observe a read with soft clipping") {
    val re = new ReadExplorer(sc.parallelize(Seq[Observation]()))

    val read = AlignmentRecord.newBuilder()
      .setStart(10L)
      .setEnd(15L)
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setMapq(40)
      .setSequence("GGGACTGAGGG")
      .setQual(":::?:::::::")
      .setCigar("3S5M3S")
      .setMismatchingPositions("0A4")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations((read, 0L))
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 5)
    assert(observations.filter(_.pos.pos != 10L).forall(_.phred == 25))
    assert(observations.forall(_.mapq == Some(40)))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.forall(_.distanceToNearestReadDeletion.isEmpty))
    assert(observations.forall(_.distanceToNearestReadInsertion.isEmpty))
    assert(observations.head.mismatchQScoreSum === Some(30))
    assert(observations.forall(_.alignedReadLen == 5))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.offsetInAlignment === 0)
    assert(observations.filter(_.pos.pos == 10L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 10L).head.phred === 30)
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.offsetInAlignment === 1)
    assert(observations.filter(_.pos.pos == 11L).head.allele === "C")
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.offsetInAlignment === 2)
    assert(observations.filter(_.pos.pos == 12L).head.allele === "T")
    assert(observations.filter(_.pos.pos == 13L).length === 1)
    assert(observations.filter(_.pos.pos == 13L).head.offsetInAlignment === 3)
    assert(observations.filter(_.pos.pos == 13L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 14L).length === 1)
    assert(observations.filter(_.pos.pos == 14L).head.offsetInAlignment === 4)
    assert(observations.filter(_.pos.pos == 14L).head.allele === "A")
  }

  sparkTest("observe a read with soft clipping and insertion") {
    val re = new ReadExplorer(sc.parallelize(Seq[Observation]()))

    val read = AlignmentRecord.newBuilder()
      .setStart(10L)
      .setEnd(15L)
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setMapq(40)
      .setSequence("GGGACTGAGG")
      .setQual(":::?:::::::")
      .setCigar("3S1M2I2M2S")
      .setMismatchingPositions("0A4")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations((read, 0L))
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 3)
    assert(observations.filter(_.pos.pos != 10L).forall(_.phred == 25))
    assert(observations.forall(_.mapq == Some(40)))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.forall(_.distanceToNearestReadDeletion.isEmpty))
    assert(observations.exists(_.distanceToNearestReadInsertion.isDefined))
    //changed
    assert(observations.head.mismatchQScoreSum === Some(25))
    assert(observations.forall(_.alignedReadLen == 5))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.distanceToNearestReadInsertion === Some(1))
    assert(observations.filter(_.pos.pos == 10L).head.allele === "ACT")
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.distanceToNearestReadInsertion === Some(-1))
    assert(observations.filter(_.pos.pos == 11L).head.allele === "G")
    // changed
    assert(observations.filter(_.pos.pos == 11L).head.phred === 25)
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.distanceToNearestReadInsertion === Some(-2))
    assert(observations.filter(_.pos.pos == 12L).head.allele === "A")
  }

  sparkTest("observe a read with hard clipping") {
    val re = new ReadExplorer(sc.parallelize(Seq[Observation]()))

    val read = AlignmentRecord.newBuilder()
      .setStart(10L)
      .setEnd(15L)
      .setContig(Contig.newBuilder()
        .setContigName("chr1")
        .build())
      .setMapq(40)
      .setSequence("ACTGA")
      .setQual("::::?")
      .setCigar("3H5M3H")
      .setMismatchingPositions("4C0")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations((read, 0L))
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 5)
    assert(observations.filter(_.pos.pos != 14L).forall(_.phred == 25))
    assert(observations.forall(_.mapq == Some(40)))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.forall(_.distanceToNearestReadDeletion.isEmpty))
    assert(observations.forall(_.distanceToNearestReadInsertion.isEmpty))
    assert(observations.head.mismatchQScoreSum === Some(30))
    assert(observations.forall(_.alignedReadLen == 5))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.offsetInAlignment === 0)
    assert(observations.filter(_.pos.pos == 10L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.offsetInAlignment === 1)
    assert(observations.filter(_.pos.pos == 11L).head.allele === "C")
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.offsetInAlignment === 2)
    assert(observations.filter(_.pos.pos == 12L).head.allele === "T")
    assert(observations.filter(_.pos.pos == 13L).length === 1)
    assert(observations.filter(_.pos.pos == 13L).head.offsetInAlignment === 3)
    assert(observations.filter(_.pos.pos == 13L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 14L).length === 1)
    assert(observations.filter(_.pos.pos == 14L).head.offsetInAlignment === 4)
    assert(observations.filter(_.pos.pos == 14L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 14L).head.phred === 30)
  }

}
