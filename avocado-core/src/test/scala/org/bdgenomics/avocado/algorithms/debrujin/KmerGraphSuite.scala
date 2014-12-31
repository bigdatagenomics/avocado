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
import scala.collection.immutable.{ NumericRange, SortedMap }
import scala.collection.mutable.ArrayBuffer

class KmerGraphSuite extends AvocadoFunSuite {

  override val properties = Map(("spark.serializer", "org.apache.spark.serializer.KryoSerializer"),
    ("spark.kryo.registrator", "org.bdgenomics.adam.serialization.ADAMKryoRegistrator"),
    ("spark.kryoserializer.buffer.mb", "128"),
    ("spark.kryo.referenceTracking", "true"))

  test("cannot build a graph without reads attached to at least one sample") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 115L)

    intercept[AssertionError] {
      val graph = KmerGraph(5, Seq((region, ref)), Seq())
    }
  }

  test("put reference into graph") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 115L)

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
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("GAGACATGC")
        .setQual("*********")
        .setRecordGroupSample("sample1")
        .setMapq(50)
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
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("GAGACATGC")
        .setQual("*********")
        .setRecordGroupSample("sample1")
        .setMapq(50)
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
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGACACATGC")
        .setQual("88888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
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
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGAGTACATGC")
        .setQual("888888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
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

  test("put reads into graph, introduce a deletion") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGAGACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGAACATGC")
        .setQual("8888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
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

  test("put reads into graph, introduce a bubble with two alleles") {
    val ref = "ACACTGAGACATGC"
    val region = ReferenceRegion("chr1", 100L, 114L)

    val graphs = KmerGraph(5, Seq((region, ref)), Seq(AlignmentRecord.newBuilder()
      .setSequence("ACACTGACACATGC")
      .setQual("88888888888888")
      .setRecordGroupSample("sample1")
      .setMapq(50)
      .build(),
      AlignmentRecord.newBuilder()
        .setSequence("ACACTGAAACATGC")
        .setQual("88888888888888")
        .setRecordGroupSample("sample1")
        .setMapq(50)
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

  /**
   * // shamelessly borrowed from the indel realigner while we are refactoring...
   * def getReferenceFromReads(reads: Seq[RichAlignmentRecord]): String = {
   * // get reference and range from a single read
   * val readRefs = reads.map((r: RichAlignmentRecord) => {
   * (r.mdTag.get.getReference(r), r.getStart.toLong to r.getEnd)
   * })
   * .toSeq
   * .sortBy(_._2.head)
   *
   * // fold over sequences and append - sequence is sorted at start
   * val ref = readRefs.reverse.foldRight[(String, Long)](("", readRefs.head._2.head))((refReads: (String, NumericRange[Long]), reference: (String, Long)) => {
   * if (refReads._2.end < reference._2) {
   * reference
   * } else if (reference._2 >= refReads._2.head) {
   * (reference._1 + refReads._1.substring((reference._2 - refReads._2.head).toInt), refReads._2.end)
   * } else {
   * // there is a gap in the sequence
   * throw new IllegalArgumentException("Current sequence has a gap at " + reference._2 + "with " + refReads._2.head + "," + refReads._2.end + ".")
   * }
   * })
   *
   * ref._1
   * }
   *
   * def na12878_chr20_snp_reads: RDD[RichAlignmentRecord] = {
   * val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058.sam").getFile
   * val reads: RDD[AlignmentRecord] = sc.loadAlignments(path)
   * reads.map(r => RichAlignmentRecord(r))
   * }
   *
   * def more_na12878_chr20_snp_reads: RDD[RichAlignmentRecord] = {
   * val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058_longer.sam").getFile
   * val reads: RDD[AlignmentRecord] = sc.loadAlignments(path)
   * reads.map(r => RichAlignmentRecord(r))
   * }
   *
   * def makeRead(sequence: String,
   * start: Long,
   * cigar: String,
   * mdtag: String,
   * length: Int,
   * qualities: Seq[Int],
   * id: Int = 0): RichAlignmentRecord = {
   *
   * val contig = Contig.newBuilder()
   * .setContigName("chr1")
   * .build()
   *
   * RichAlignmentRecord(AlignmentRecord.newBuilder()
   * .setReadName("read" + id.toString)
   * .setStart(start)
   * .setReadMapped(true)
   * .setCigar(cigar)
   * .setSequence(sequence)
   * .setReadNegativeStrand(false)
   * .setMapq(60)
   * .setQual(qualities.map(_.toChar.toString).reduce(_ + _))
   * .setMismatchingPositions(mdtag)
   * .setContig(contig)
   * .setRecordGroupSample("sample")
   * .build())
   * }
   *
   * test("Test the creation of several haplotype strings.") {
   * val kms = ArrayBuffer[Kmer](new Kmer("ACG"),
   * Kmer("CGA"),
   * Kmer("GAG"))
   * val kp = new KmerPath(kms)
   *
   * assert(kp.haplotypeString === "ACGAG")
   * }
   *
   * test("put kmers into graph for a single, simple read") {
   * val reference = "TACCAAT"
   * val read = makeRead("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
   *
   * val graph = KmerGraph(3, 7, 7, reference, 3, 5)
   * graph.insertRead(read)
   * assert(graph.kmers.size === 5)
   * assert(graph.prefixSet.contains("TA"))
   * assert(graph.prefixSet.contains("AC"))
   * assert(graph.prefixSet.contains("CC"))
   * assert(graph.prefixSet.contains("CA"))
   *
   * assert(graph.kmerSequences.toSet.contains("TAC"))
   * assert(graph.kmerSequences.toSet.contains("ACC"))
   * assert(graph.kmerSequences.toSet.contains("CCA"))
   * assert(graph.kmerSequences.toSet.contains("CAA"))
   *
   * assert(graph.allPaths.size === 1)
   * assert(graph.allPaths.head.haplotypeString === "TACCAAT")
   * assert(graph.allPaths.head.weight === (reference.length - 2))
   * }
   *
   * test("put kmers into graph for a small set of reads without a polymorphism") {
   * val reference = "TACCAATGTAA"
   * val read0 = makeRead("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
   * val read1 = makeRead("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
   * val read2 = makeRead("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
   * val read3 = makeRead("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
   * val read4 = makeRead("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
   *
   * val readBucket = Seq(read0, read1, read2, read3, read4)
   * val kmerGraph = KmerGraph(4, 7, 7, reference, readBucket, 4)
   *
   * assert(kmerGraph.sourceKmer.prefix === "TAC")
   * assert(kmerGraph.sourceKmer.suffix === 'C')
   * assert(kmerGraph.sinkKmer.prefix === "GTA")
   * assert(kmerGraph.sinkKmer.suffix === 'A')
   * /*
   * Expect following Kmers:
   * - TAC[C]
   * - ACC[A]
   * - CCA[A]
   * - CAA[T]
   * - AAT[G]
   * - ATG[T]
   * - TGT[A]
   * - GTA[A]
   * */
   * assert(kmerGraph.prefixSet.size === 8)
   * assert(kmerGraph.kmers.size === 8)
   * assert(kmerGraph.prefixSet.contains("TAC"))
   * assert(kmerGraph.prefixSet.contains("ACC"))
   * assert(kmerGraph.prefixSet.contains("CCA"))
   * assert(kmerGraph.prefixSet.contains("CAA"))
   * assert(kmerGraph.prefixSet.contains("AAT"))
   * assert(kmerGraph.prefixSet.contains("ATG"))
   * assert(kmerGraph.prefixSet.contains("TGT"))
   * assert(kmerGraph.prefixSet.contains("GTA"))
   *
   * assert(kmerGraph.allPaths.size === 1)
   * assert(kmerGraph.allPaths.head.haplotypeString === "TACCAATGTAA")
   * // middle kmers are covered 4x, dropping off to 1x at ends
   * assert(kmerGraph.allPaths.head.weight === (2 * 4 + 2 * 3 + 2 * 2 + 2 * 1))
   * }
   *
   * test("merge simple kmer sequence") {
   * val seqs1 = KmerGraph.buildPrefixMap(Seq("AAC", "ACT", "CTG"))
   * val seqs2 = KmerGraph.buildPrefixMap(Seq("TAC", "ACT", "CTC"))
   *
   * assert(seqs1.flatMap(kv => kv._2).size === 3)
   * assert(seqs2.flatMap(kv => kv._2).size === 3)
   *
   * val merge = KmerGraph.mergeGraphs(seqs1, seqs2)
   *
   * assert(merge.flatMap(kv => kv._2).size === 5)
   * }
   *
   * test("put kmers into graph for a small set of reads with a polymorphism") {
   * val reference = "TACCAATGTAA"
   * val read0 = makeRead("TACCCAT", 0L, "7M", "4A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
   * val read1 = makeRead("ACCCATG", 1L, "7M", "3A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
   * val read2 = makeRead("CCCATGT", 2L, "7M", "2A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
   * val read3 = makeRead("CCATGTA", 3L, "7M", "1A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
   * val read4 = makeRead("CATGTAA", 4L, "7M", "0A6", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
   * val read5 = makeRead("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
   * val read6 = makeRead("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
   * val read7 = makeRead("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
   * val read8 = makeRead("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
   * val read9 = makeRead("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)
   *
   * val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
   * val kmerGraph = KmerGraph(4, 7, 7, reference, readBucket, 4)
   *
   * /*
   * Expect following Kmers:
   * - TAC[C]
   * - ACC[C,A]
   * - CCC[A]
   * - CCA[T,A]
   * - CAA[T]
   * - CAT[G]
   * - AAT[G]
   * - ATG[T]
   * - TGT[A]
   * - GTA[A]
   * */
   * assert(kmerGraph.prefixSet.size === 10)
   * assert(kmerGraph.kmers.size === 12)
   * assert(kmerGraph.prefixSet.contains("TAC"))
   * assert(kmerGraph.prefixSet.contains("ACC"))
   * assert(kmerGraph.prefixSet.contains("CCC"))
   * assert(kmerGraph.prefixSet.contains("CCA"))
   * assert(kmerGraph.prefixSet.contains("CAA"))
   * assert(kmerGraph.prefixSet.contains("CAT"))
   * assert(kmerGraph.prefixSet.contains("AAT"))
   * assert(kmerGraph.prefixSet.contains("ATG"))
   * assert(kmerGraph.prefixSet.contains("TGT"))
   * assert(kmerGraph.prefixSet.contains("GTA"))
   *
   * assert(kmerGraph.allPaths.size === 4)
   * val haplotypes = kmerGraph.allPaths.map(_.haplotypeString)
   * assert(haplotypes.contains("TACCAATGTAA"))
   * assert(haplotypes.contains("TACCCATGTAA"))
   * assert(haplotypes.contains("TACCCAATGTAA"))
   * assert(haplotypes.contains("TACCATGTAA"))
   * }
   *
   * sparkTest("put reads into graph for real data") {
   * val reads = na12878_chr20_snp_reads.collect.toSeq
   *
   * val reference = getReferenceFromReads(reads)
   *
   * val kmerGraph = KmerGraph(20, 101, reference.length, reference, reads, 40)
   * assert(kmerGraph.allPaths.size === 30)
   * assert(kmerGraph.kmers.size === 678)
   *
   * kmerGraph.removeSpurs()
   * assert(kmerGraph.allPaths.size === 30)
   * assert(kmerGraph.kmers.size === 409)
   * }
   *
   * test("check spur trimming function on simple graph") {
   *
   * val reference = "TACCAATGTAA"
   * val read0 = makeRead("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
   * val read1 = makeRead("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
   * val read2 = makeRead("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
   * val read3 = makeRead("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
   * val read4 = makeRead("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
   * val read5 = makeRead("CCAATAT", 2L, "7M", "5G1", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
   *
   * val readBucket = Seq(read0, read1, read2, read3, read4, read5)
   * val kmerGraph = KmerGraph(4, 7, 7, reference, readBucket, 4)
   *
   * assert(kmerGraph.allPaths.size === 1)
   * assert(kmerGraph.kmers.size === 10)
   *
   * kmerGraph.removeSpurs()
   *
   * assert(kmerGraph.allPaths.size === 1)
   * assert(kmerGraph.kmers.size === 8)
   * }
   *
   * sparkTest("put reads into graph for larger real dataset") {
   * val reads = more_na12878_chr20_snp_reads.collect.toSeq
   *
   * val reference = getReferenceFromReads(reads)
   *
   * val kmerGraph = KmerGraph(20, 101, reference.length, reference, reads, 40)
   *
   * assert(kmerGraph.kmers.size === 3271)
   *
   * kmerGraph.removeSpurs()
   *
   * assert(kmerGraph.kmers.size === 1714)
   *
   * kmerGraph.trimLowCoverageKmers(0.05)
   * kmerGraph.removeSpurs()
   *
   * assert(kmerGraph.kmers.size === 659)
   * assert(kmerGraph.allPaths.size === 12)
   * }
   *
   * test("generate a limited number of haplotypes for a graph with repeats") {
   * val reference = "TACCAAAAATGTAA"
   * val read0 = makeRead("TACCAAAAATGTAA", 0L, "15M", "15", 15,
   * Seq(50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50), 0)
   * val kmerGraph = KmerGraph(4, 15, 15, reference, Seq(read0), 4)
   * assert(kmerGraph.allPaths.size === 7)
   * }
   */
}
