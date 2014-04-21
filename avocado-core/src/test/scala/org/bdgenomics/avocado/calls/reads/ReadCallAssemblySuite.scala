/*
 * Copyright (c) 2013-2014. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

package org.bdgenomics.avocado.calls.reads

import org.bdgenomics.adam.avro.{ ADAMRecord, ADAMGenotypeAllele }
import org.bdgenomics.adam.models.ADAMVariantContext
import org.bdgenomics.avocado.algorithms.debrujin._
import org.bdgenomics.avocado.algorithms.hmm._
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.scalatest.FunSuite
import scala.collection.JavaConversions._
import scala.collection.mutable.{ ArrayBuffer, Buffer, HashMap, HashSet, PriorityQueue, StringBuilder }

class ReadCallAssemblySuite extends FunSuite {

  val rcap_short = new ReadCallAssemblyPhaser(4, 0)

  def make_read(sequence: String,
                start: Long,
                cigar: String,
                mdtag: String,
                length: Int,
                qualities: Seq[Int],
                id: Int = 0): ADAMRecord = {
    ADAMRecord.newBuilder()
      .setReadName("read" + id.toString)
      .setStart(start)
      .setReadMapped(true)
      .setCigar(cigar)
      .setSequence(sequence)
      .setReadNegativeStrand(false)
      .setMapq(60)
      .setQual(qualities.map(_.toChar.toString).reduce(_ + _))
      .setMismatchingPositions(mdtag)
      .setReferenceName("chr1")
      .setReferenceId(1)
      .setRecordGroupSample("sample")
      .build()
  }

  test("Test the creation of several haplotype strings.") {
    val kms = ArrayBuffer[Kmer](new Kmer(new KmerPrefix("AC"), 'G'),
      new Kmer(new KmerPrefix("CG"), 'A'),
      new Kmer(new KmerPrefix("GA"), 'G'))
    val kp = new KmerPath(kms)

    assert(kp.asHaplotypeString === "ACGAG")
  }

  test("Test log sum for similar values") {
    val hp = new HaplotypePair(null, null)
    val sum = hp.exactLogSumExp10(1.0, 1.0)

    assert(1.3 * 0.99 < sum && 1.3 * 1.01 > sum)
  }

  test("Test log sum for dissimilar values") {
    val hp = new HaplotypePair(null, null)
    val sum = hp.exactLogSumExp10(1.0, -3.0)

    assert(1.00004342 * 0.99 < sum && 1.00004342 * 1.01 > sum)
  }

  test("put kmers into graph for a single, simple read") {
    val reference = "TACCAAT"
    val read = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val aread = AssemblyRead(read)

    val graph = rcap_short.assemble(Seq(read), reference)

    assert(graph.kmers.size === 4)
    assert(graph.prefixes.contains("TAC"))
    assert(graph.prefixes.contains("ACC"))
    assert(graph.prefixes.contains("CCA"))
    assert(graph.prefixes.contains("CAA"))
    assert(graph.kmers.values.toList.flatMap(k => k).size == 4)

    assert(graph.allPaths.length === 1)
    assert(graph.allPaths.head.asHaplotypeString === "TACCAAT")
    assert(graph.allPaths.head.multSum === (reference.length - 3))
  }

  test("put kmers into graph for a small set of reads without a polymorphism") {
    val reference = "TACCAATGTAA"
    val read0 = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)

    val readBucket = Seq(read0, read1, read2, read3, read4)
    val kmerGraph = rcap_short.assemble(readBucket, reference)

    assert(kmerGraph.sourceKmer.prefix.string === "TAC")
    assert(kmerGraph.sourceKmer.suffix === 'C')
    assert(kmerGraph.sinkKmer.prefix.string === "GTA")
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
    assert(kmerGraph.prefixes.size === 8)
    assert(kmerGraph.kmers.values.toList.flatMap(v => v).size === 8)
    assert(kmerGraph.prefixes.contains("TAC"))
    assert(kmerGraph.prefixes.contains("ACC"))
    assert(kmerGraph.prefixes.contains("CCA"))
    assert(kmerGraph.prefixes.contains("CAA"))
    assert(kmerGraph.prefixes.contains("AAT"))
    assert(kmerGraph.prefixes.contains("ATG"))
    assert(kmerGraph.prefixes.contains("TGT"))
    assert(kmerGraph.prefixes.contains("GTA"))

    assert(kmerGraph.allPaths.length === 1)
    assert(kmerGraph.allPaths.head.asHaplotypeString === "TACCAATGTAA")
    // middle kmers are covered 4x, dropping off to 1x at ends
    assert(kmerGraph.allPaths.head.multSum === (2 * 4 + 2 * 3 + 2 * 2 + 2 * 1))
  }

  test("put kmers into graph for a small set of reads with a polymorphism") {
    val reference = "TACCAATGTAA"
    val read0 = make_read("TACCCAT", 0L, "7M", "4A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCCATG", 1L, "7M", "3A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCCATGT", 2L, "7M", "2A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CCATGTA", 3L, "7M", "1A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("CATGTAA", 4L, "7M", "0A6", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
    val read5 = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
    val read6 = make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
    val read7 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
    val read8 = make_read("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
    val read9 = make_read("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
    val kmerGraph = rcap_short.assemble(readBucket, reference)

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
    assert(kmerGraph.prefixes.size === 10)
    assert(kmerGraph.kmers.values.toList.flatMap(v => v).size === 12)
    assert(kmerGraph.prefixes.contains("TAC"))
    assert(kmerGraph.prefixes.contains("ACC"))
    assert(kmerGraph.prefixes.contains("CCC"))
    assert(kmerGraph.prefixes.contains("CCA"))
    assert(kmerGraph.prefixes.contains("CAA"))
    assert(kmerGraph.prefixes.contains("CAT"))
    assert(kmerGraph.prefixes.contains("AAT"))
    assert(kmerGraph.prefixes.contains("ATG"))
    assert(kmerGraph.prefixes.contains("TGT"))
    assert(kmerGraph.prefixes.contains("GTA"))

    assert(kmerGraph.allPaths.length === 4)
    val haplotypes = kmerGraph.allPaths.map(_.asHaplotypeString)
    assert(haplotypes.contains("TACCAATGTAA"))
    assert(haplotypes.contains("TACCCATGTAA"))
    assert(haplotypes.contains("TACCCAATGTAA"))
    assert(haplotypes.contains("TACCATGTAA"))
  }

  test("\"Call\" hom ref, ~10x coverage") {
    val reference = "TACCAATGTAA"
    val read0 = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
    val read5 = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
    val read6 = make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
    val read7 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
    val read8 = make_read("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
    val read9 = make_read("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
    val kmerGraph = rcap_short.assemble(readBucket, reference)

    val variants = rcap_short.phaseAssembly(readBucket, kmerGraph, reference)

    assert(variants.length === 0)
  }

  test("Call simple het SNP, ~10x coverage") {
    val reference = "TACCAATGTAA"
    val read0 = make_read("TACCCAT", 0L, "7M", "4A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCCATG", 1L, "7M", "3A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCCATGT", 2L, "7M", "2A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CCATGTA", 3L, "7M", "1A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("CATGTAA", 4L, "7M", "0A6", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
    val read5 = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
    val read6 = make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
    val read7 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
    val read8 = make_read("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
    val read9 = make_read("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
    val kmerGraph = rcap_short.assemble(readBucket, reference)

    val variants = rcap_short.phaseAssembly(readBucket, kmerGraph, reference)

    assert(variants.length === 1)
    assert(variants.head.position.pos === 4L)
    assert(variants.head.variant.variant.getReferenceAllele === "A")
    assert(variants.head.variant.variant.getVariantAllele === "C")
    val alleles: List[ADAMGenotypeAllele] = asScalaBuffer(variants.head.genotypes.head.getAlleles).toList
    assert(alleles.length === 2)
    assert(alleles.head === ADAMGenotypeAllele.Ref)
    assert(alleles.last === ADAMGenotypeAllele.Alt)
  }

  test("Call simple hom SNP, ~10x coverage") {
    val reference = "TACCAATGTAA"
    val read0 = make_read("TACCCAT", 0L, "7M", "4A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCCATG", 1L, "7M", "3A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCCATGT", 2L, "7M", "2A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CCATGTA", 3L, "7M", "1A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("CATGTAA", 4L, "7M", "0A6", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
    val read5 = make_read("TACCCAT", 0L, "7M", "4A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
    val read6 = make_read("ACCCATG", 1L, "7M", "3A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
    val read7 = make_read("CCCATGT", 2L, "7M", "2A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
    val read8 = make_read("CCATGTA", 3L, "7M", "1A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
    val read9 = make_read("CATGTAA", 4L, "7M", "0A6", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
    val kmerGraph = rcap_short.assemble(readBucket, reference)

    val variants = rcap_short.phaseAssembly(readBucket, kmerGraph, reference)

    assert(variants.length === 1)
    assert(variants.head.position.pos === 4L)
    assert(variants.head.variant.variant.getReferenceAllele === "A")
    assert(variants.head.variant.variant.getVariantAllele === "C")
    val alleles: List[ADAMGenotypeAllele] = asScalaBuffer(variants.head.genotypes.head.getAlleles).toList
    assert(alleles.length === 2)
    assert(alleles.head === ADAMGenotypeAllele.Alt)
    assert(alleles.last === ADAMGenotypeAllele.Alt)
  }

  test("Call simple het INS, ~10x coverage") {
    val reference = "TACCAATGTAA"
    val read0 = make_read("TACCCAA", 0L, "4M1I2M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCCAAT", 1L, "3M1I3M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCCAATG", 2L, "2M1I4M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CCCAATG", 2L, "2M1I4M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("CCAATGT", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
    val read5 = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
    val read6 = make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
    val read7 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
    val read8 = make_read("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
    val read9 = make_read("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
    val kmerGraph = rcap_short.assemble(readBucket, reference)

    val variants = rcap_short.phaseAssembly(readBucket, kmerGraph, reference)

    assert(variants.length === 1)
    assert(variants.head.position.pos === 2L)
    assert(variants.head.variant.variant.getReferenceAllele === "C")
    assert(variants.head.variant.variant.getVariantAllele === "CC")
    val alleles: List[ADAMGenotypeAllele] = asScalaBuffer(variants.head.genotypes.head.getAlleles).toList
    assert(alleles.length === 2)
    assert(alleles.head === ADAMGenotypeAllele.Ref)
    assert(alleles.last === ADAMGenotypeAllele.Alt)
  }

  test("Call simple het DEL, ~10x coverage") {
    val reference = "TACCAATGTAA"
    val read0 = make_read("TACCATG", 0L, "4M1D2M", "4^A2", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCATGT", 1L, "3M1D3M", "3^A3", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCATGTA", 2L, "2M1D4M", "2^A4", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CATGTAA", 3L, "1M1D5M", "1^A5", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("TACCAAT", 0L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
    val read5 = make_read("ACCAATG", 1L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 5)
    val read6 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 6)
    val read7 = make_read("CCAATGT", 2L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 7)
    val read8 = make_read("CAATGTA", 3L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 8)
    val read9 = make_read("AATGTAA", 4L, "7M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 9)

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8, read9)
    val kmerGraph = rcap_short.assemble(readBucket, reference)

    val variants = rcap_short.phaseAssembly(readBucket, kmerGraph, reference)

    assert(variants.length === 1)
    assert(variants.head.position.pos === 4L)
    assert(variants.head.variant.variant.getReferenceAllele === "CA")
    assert(variants.head.variant.variant.getVariantAllele === "C")
    val alleles: List[ADAMGenotypeAllele] = asScalaBuffer(variants.head.genotypes.head.getAlleles).toList
    assert(alleles.length === 2)
    assert(alleles.head === ADAMGenotypeAllele.Ref)
    assert(alleles.last === ADAMGenotypeAllele.Alt)
  }
}
