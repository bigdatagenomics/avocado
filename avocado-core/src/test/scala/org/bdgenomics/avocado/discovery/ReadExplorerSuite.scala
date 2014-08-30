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

import org.bdgenomics.adam.util.SparkFunSuite
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.formats.avro.{ AlignmentRecord, Contig }

class ReadExplorerSuite extends SparkFunSuite {

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
      .setQual(":::::")
      .setCigar("5M")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations(read)
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 5)
    assert(observations.forall(_.phred == 25))
    assert(observations.forall(_.mapq == 40))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
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
      .setQual(":::::")
      .setCigar("2M2D3M")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations(read)
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 7)
    assert(observations.filter(_.allele != "_").forall(_.phred == 25))
    assert(observations.filter(_.allele == "_").forall(_.phred == 40))
    assert(observations.forall(_.mapq == 40))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.allele === "A")
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.allele === "C")
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.allele === "_")
    assert(observations.filter(_.pos.pos == 13L).length === 1)
    assert(observations.filter(_.pos.pos == 13L).head.allele === "_")
    assert(observations.filter(_.pos.pos == 14L).length === 1)
    assert(observations.filter(_.pos.pos == 14L).head.allele === "T")
    assert(observations.filter(_.pos.pos == 15L).length === 1)
    assert(observations.filter(_.pos.pos == 15L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 16L).length === 1)
    assert(observations.filter(_.pos.pos == 16L).head.allele === "A")
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
      .setQual(":::::")
      .setCigar("1M2I2M")
      .setRecordGroupSample("sample1")
      .build()

    val observations = re.readToObservations(read)
      .flatMap(o => o match {
        case ao: AlleleObservation => Some(ao)
        case _                     => None
      })

    assert(observations.length === 3)
    assert(observations.forall(_.phred == 25))
    assert(observations.forall(_.mapq == 40))
    assert(observations.forall(_.pos.referenceName == "chr1"))
    assert(observations.forall(_.sample == "sample1"))
    assert(observations.forall(!_.onNegativeStrand))
    assert(observations.filter(_.pos.pos == 10L).length === 1)
    assert(observations.filter(_.pos.pos == 10L).head.allele === "ACT")
    assert(observations.filter(_.pos.pos == 11L).length === 1)
    assert(observations.filter(_.pos.pos == 11L).head.allele === "G")
    assert(observations.filter(_.pos.pos == 12L).length === 1)
    assert(observations.filter(_.pos.pos == 12L).head.allele === "A")
  }
}
