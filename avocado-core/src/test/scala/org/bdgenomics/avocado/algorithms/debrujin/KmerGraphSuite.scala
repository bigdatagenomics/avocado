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
package org.bdgenomics.avocado.algorithms.debrujin

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichAlignmentRecord
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.avocado.models.AlleleObservation
import org.bdgenomics.formats.avro.{ Contig, AlignmentRecord }
import org.scalatest.FunSuite
import scala.collection.immutable.{ NumericRange, SortedMap }
import scala.collection.mutable.ArrayBuffer

class KmerGraphSuite extends AvocadoFunSuite {

  test("cannot build a graph without reads attached to at least one sample") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 115L)

    intercept[AssertionError] {
      val graph = KmerGraph(5, Seq((region, ref)), Seq())
    }
  }

  test("put reference into graph") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("")
      .setQual("")
      .setRecordGroupSample("sample1")
      .setMapq(0)
      .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 10)
    assert(graph.nonRefSize === 0)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 0)

    val observations = graph.toObservations

    assert(observations.size === 10)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
  }

  test("put reads into graph, all match reference") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGA")
      .setQual("*********")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("GAGACATGC")
        .setQual("*********")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 10)
    assert(graph.nonRefSize === 0)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 0)

    val observations = graph.toObservations

    assert(observations.size === 20)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 10)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
  }

  test("put reads into graph, contains a spur") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGT")
      .setQual("*********")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("GAGACATGC")
        .setQual("*********")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 11)
    assert(graph.nonRefSize === 1)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 1)

    val observations = graph.toObservations

    assert(observations.size === 19)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 9)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
  }

  test("put reads into graph, introduce a bubble") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGACACATGC")
        .setQual("88888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 15)
    assert(graph.nonRefSize === 5)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 0)

    val observations = graph.toObservations

    assert(observations.size === 30)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 20)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos != 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos == 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 2))
  }

  test("put reads into graph, introduce a bubble with an insert") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGAGTACATGC")
        .setQual("888888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 15)
    assert(graph.nonRefSize === 5)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 0)

    val observations = graph.toObservations

    assert(observations.size === 31)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 21)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos != 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos == 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 2))
  }

  test("put reads into graph, introduce a bubble with a complex insert") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGATTAACATGC")
        .setQual("8888888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 17)
    assert(graph.nonRefSize === 7)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 0)

    val observations = graph.toObservations

    assert(observations.size === 32)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 22)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos != 107L && kv._1.pos != 106L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos == 107L || kv._1.pos == 106L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 2))
  }

  test("put reads into graph, introduce a deletion") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGAACATGC")
        .setQual("8888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 14)
    assert(graph.nonRefSize === 4)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 0)

    val observations = graph.toObservations

    assert(observations.size === 29)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 19)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos != 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos == 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 2))
  }

  test("put reads into graph, introduce a complex deletion") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGCCATGC")
        .setQual("888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 15)
    assert(graph.nonRefSize === 5)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)
    assert(graph.spurs === 0)

    val observations = graph.toObservations

    assert(observations.size === 30)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 20)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos != 105L && kv._1.pos != 106L && kv._1.pos != 108L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos == 105L || kv._1.pos == 106L || kv._1.pos == 108L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 2))
  }

  test("put reads into graph, introduce a bubble with two alleles") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGACACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGAAACATGC")
        .setQual("88888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 20)
    assert(graph.nonRefSize === 10)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)

    val observations = graph.toObservations

    assert(observations.size === 30)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 20)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos != 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
    observations.groupBy(_.pos)
      .filter(kv => kv._1.pos == 107L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 3))
  }

  test("put reads into graph, introduce a mnp bubble") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .setReadNegativeStrand(false)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTCAAACATGC")
        .setQual("88888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
        .setReadNegativeStrand(false)
        .build()))

    assert(graphs.size === 1)
    val graph = graphs.head

    assert(graph.sample === "sample1")
    assert(graph.size === 17)
    assert(graph.nonRefSize === 7)
    assert(graph.sources === 1)
    assert(graph.sinks === 1)

    val observations = graph.toObservations

    assert(observations.size === 30)
    observations.foreach(o => {
      val r = o.pos
      assert(r.referenceName === "chr1")
      assert(r.pos >= 100L && r.pos < 115L)
    })
    assert(observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).size === 20)
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .foreach(kv => assert(kv._2.length === 2))
    observations.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    }).groupBy(_.pos)
      .filter(kv => kv._1.pos != 107L && kv._1.pos != 105L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 1))
    observations.groupBy(_.pos)
      .filter(kv => kv._1.pos == 107L || kv._1.pos == 105L)
      .map(kv => kv._2
        .map(ao => ao.allele)
        .toSet).foreach(v => assert(v.size === 2))
  }
}
