/*
 * Copyright (c) 2014. Mount Sinai School of Medicine
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use file except in compliance with the License.
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

package org.bdgenomics.avocado.algorithms.debrujin

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.avro.{ ADAMContig, ADAMRecord }
import org.bdgenomics.adam.rdd.ADAMContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichADAMRecord
import org.bdgenomics.adam.util.SparkFunSuite
import org.bdgenomics.avocado.calls.reads.ReadCallAssemblyPhaser
import parquet.filter.UnboundRecordFilter
import scala.collection.mutable.ArrayBuffer

class KmerGraphSuite extends SparkFunSuite {

  def na12878_chr20_snp_reads: RDD[RichADAMRecord] = {
    val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058.sam").getFile
    sc.adamLoad[ADAMRecord, UnboundRecordFilter](path).map(r => RichADAMRecord(r))
  }

  def makeRead(sequence: String,
               start: Long,
               cigar: String,
               mdtag: String,
               length: Int,
               qualities: Seq[Int],
               id: Int = 0): RichADAMRecord = {

    val contig = ADAMContig.newBuilder()
      .setContigName("chr1")
      .build()

    RichADAMRecord(ADAMRecord.newBuilder()
      .setReadName("read" + id.toString)
      .setStart(start)
      .setReadMapped(true)
      .setCigar(cigar)
      .setSequence(sequence)
      .setReadNegativeStrand(false)
      .setMapq(60)
      .setQual(qualities.map(_.toChar.toString).reduce(_ + _))
      .setMismatchingPositions(mdtag)
      .setContig(contig)
      .setRecordGroupSample("sample")
      .build())
  }

  test("Test the creation of several haplotype strings.") {
    val kms = ArrayBuffer[Kmer](new Kmer("ACG"),
      Kmer("CGA"),
      Kmer("GAG"))
    val kp = new KmerPath(kms)

    assert(kp.haplotypeString === "ACGAG")
  }

  test("put kmers into graph for a single, simple read") {
    val reference = "TACCAAT"
    val read = makeRead("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)

    val graph = KmerGraph(3, 7, 7, reference, 3)
    graph.insertRead(read)
    assert(graph.kmers.size === 5)
    assert(graph.prefixSet.contains("TA"))
    assert(graph.prefixSet.contains("AC"))
    assert(graph.prefixSet.contains("CC"))
    assert(graph.prefixSet.contains("CA"))

    assert(graph.kmerSequences.toSet.contains("TAC"))
    assert(graph.kmerSequences.toSet.contains("ACC"))
    assert(graph.kmerSequences.toSet.contains("CCA"))
    assert(graph.kmerSequences.toSet.contains("CAA"))

    assert(graph.allPaths.size === 1)
    assert(graph.allPaths.head.haplotypeString === "TACCAAT")
    assert(graph.allPaths.head.weight === (reference.length - 2))
  }

  test("put kmers into graph for a small set of reads without a polymorphism") {
    val reference = "TACCAATGTAA"
    val read0 = makeRead("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = makeRead("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = makeRead("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = makeRead("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = makeRead("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)

    val readBucket = Seq(read0, read1, read2, read3, read4)
    val kmerGraph = KmerGraph(4, 7, 7, reference, readBucket, 4)

    assert(kmerGraph.sourceKmer.prefix === "TAC")
    assert(kmerGraph.sourceKmer.suffix === 'C')
    assert(kmerGraph.sinkKmer.prefix === "GTA")
    assert(kmerGraph.sinkKmer.suffix === 'A')
    /*
     * Expect following Kmers:
     * - TAC[C]
     * - ACC[A]
     * - CCA[A]
     * - CAA[T]
     * - AAT[G]
     * - ATG[T]
     * - TGT[A]
     * - GTA[A]
     */
    assert(kmerGraph.prefixSet.size === 8)
    assert(kmerGraph.kmers.size === 8)
    assert(kmerGraph.prefixSet.contains("TAC"))
    assert(kmerGraph.prefixSet.contains("ACC"))
    assert(kmerGraph.prefixSet.contains("CCA"))
    assert(kmerGraph.prefixSet.contains("CAA"))
    assert(kmerGraph.prefixSet.contains("AAT"))
    assert(kmerGraph.prefixSet.contains("ATG"))
    assert(kmerGraph.prefixSet.contains("TGT"))
    assert(kmerGraph.prefixSet.contains("GTA"))

    assert(kmerGraph.allPaths.size === 1)
    assert(kmerGraph.allPaths.head.haplotypeString === "TACCAATGTAA")
    // middle kmers are covered 4x, dropping off to 1x at ends
    assert(kmerGraph.allPaths.head.weight === (2 * 4 + 2 * 3 + 2 * 2 + 2 * 1))
  }

  test("put kmers into graph for a small set of reads with a polymorphism") {
    val reference = "TACCAATGTAA"
    val read0 = makeRead("TACCCAT", 0L, "7M", "4A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = makeRead("ACCCATG", 1L, "7M", "3A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = makeRead("CCCATGT", 2L, "7M", "2A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = makeRead("CCATGTA", 3L, "7M", "1A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = makeRead("CATGTAA", 4L, "7M", "0A6", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
    val read5 = makeRead("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
    val read6 = makeRead("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
    val read7 = makeRead("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
    val read8 = makeRead("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
    val read9 = makeRead("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
    val kmerGraph = KmerGraph(4, 7, 7, reference, readBucket, 4)

    /*
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
       */
    assert(kmerGraph.prefixSet.size === 10)
    assert(kmerGraph.kmers.size === 12)
    assert(kmerGraph.prefixSet.contains("TAC"))
    assert(kmerGraph.prefixSet.contains("ACC"))
    assert(kmerGraph.prefixSet.contains("CCC"))
    assert(kmerGraph.prefixSet.contains("CCA"))
    assert(kmerGraph.prefixSet.contains("CAA"))
    assert(kmerGraph.prefixSet.contains("CAT"))
    assert(kmerGraph.prefixSet.contains("AAT"))
    assert(kmerGraph.prefixSet.contains("ATG"))
    assert(kmerGraph.prefixSet.contains("TGT"))
    assert(kmerGraph.prefixSet.contains("GTA"))

    assert(kmerGraph.allPaths.size === 4)
    val haplotypes = kmerGraph.allPaths.map(_.haplotypeString)
    assert(haplotypes.contains("TACCAATGTAA"))
    assert(haplotypes.contains("TACCCATGTAA"))
    assert(haplotypes.contains("TACCCAATGTAA"))
    assert(haplotypes.contains("TACCATGTAA"))
  }

  sparkTest("put reads into graph for real data") {
    val reads = na12878_chr20_snp_reads.collect.toSeq

    // make generic assembler to get reference recovery method
    val rcap = new ReadCallAssemblyPhaser(4, 0, 4)
    val reference = rcap.getReference(reads)

    val kmerGraph = KmerGraph(20, 101, reference.length, reference, reads, 40)
    assert(kmerGraph.allPaths.size === 30)
  }
}
