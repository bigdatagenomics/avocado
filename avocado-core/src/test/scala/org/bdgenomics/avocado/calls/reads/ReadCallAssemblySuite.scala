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

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichADAMRecord
import org.bdgenomics.adam.util.SparkFunSuite
import org.bdgenomics.adam.avro.{ ADAMGenotypeAllele, ADAMRecord, ADAMContig }
import org.bdgenomics.avocado.algorithms.hmm._
import org.bdgenomics.avocado.partitioners.PartitionSet
import parquet.filter.UnboundRecordFilter
import scala.collection.JavaConversions._
import scala.collection.immutable.SortedMap

class ReadCallAssemblySuite extends SparkFunSuite {

  val emptyPartition = new PartitionSet(SortedMap[ReferenceRegion, Int]())
  val rcap_short = new ReadCallAssemblyPhaser(emptyPartition, 4, 4)
  val rcap_long = new ReadCallAssemblyPhaser(emptyPartition, 20, 40)
  val rcap_long_trim = new ReadCallAssemblyPhaser(emptyPartition, 20, 40, 5, true, Some(0.05))

  def na12878_chr20_snp_reads: RDD[RichADAMRecord] = {
    val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058.sam").getFile
    val rdd: RDD[ADAMRecord] = sc.adamLoad(path)
    rdd.map(r => RichADAMRecord(r))
  }

  def make_read(sequence: String,
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
      .setRecordGroupSample("sample")
      .setContig(contig)
      .build())
  }

  test("Test log sum for similar values") {
    val sum = HaplotypePair.exactLogSumExp10(1.0, 1.0)
    assert(1.3 * 0.99 < sum && 1.3 * 1.01 > sum)
  }

  test("Test log sum for dissimilar values") {
    val sum = HaplotypePair.exactLogSumExp10(1.0, -3.0)
    assert(1.00004342 * 0.99 < sum && 1.00004342 * 1.01 > sum)
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

    val readBucket = Seq(read0, read1, read2, read3, read4, read5, read6, read7, read8)
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
    val read0 = make_read("TACCAAA", 0L, "4M1I2M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 0)
    val read1 = make_read("ACCAAAT", 1L, "3M1I3M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 1)
    val read2 = make_read("CCAAATG", 2L, "2M1I4M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 2)
    val read3 = make_read("CAAATGT", 2L, "1M1I5M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 3)
    val read4 = make_read("AAATGTA", 3L, "1I6M", "7", 7, Seq(50, 50, 50, 50, 50, 50, 50), 4)
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
    assert(variants.head.variant.variant.getVariantAllele === "AA")
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

  sparkTest("call A->G snp on NA12878 chr20 @ 225058") {
    val reads = na12878_chr20_snp_reads.collect.toSeq
    val reference = rcap_long.getReference(reads)

    val kmerGraph = rcap_long.assemble(reads, reference)

    val variants = rcap_long.phaseAssembly(reads, kmerGraph, reference)

    assert(variants.length === 1)
    assert(variants.head.position.pos === 225057L)
    assert(variants.head.variant.variant.getReferenceAllele === "A")
    assert(variants.head.variant.variant.getVariantAllele === "G")
    val alleles: List[ADAMGenotypeAllele] = asScalaBuffer(variants.head.genotypes.head.getAlleles).toList
    assert(alleles.length === 2)
    assert(alleles.head === ADAMGenotypeAllele.Ref)
    assert(alleles.last === ADAMGenotypeAllele.Alt)
  }

  sparkTest("call doesn't change after \"aggressive\" trimming") {
    val reads = na12878_chr20_snp_reads.collect.toSeq
    val reference = rcap_long_trim.getReference(reads)

    val kmerGraph = rcap_long_trim.assemble(reads, reference)

    val variants = rcap_long_trim.phaseAssembly(reads, kmerGraph, reference)

    assert(variants.length === 1)
    assert(variants.head.position.pos === 225057L)
    assert(variants.head.variant.variant.getReferenceAllele === "A")
    assert(variants.head.variant.variant.getVariantAllele === "G")
    val alleles: List[ADAMGenotypeAllele] = asScalaBuffer(variants.head.genotypes.head.getAlleles).toList
    assert(alleles.length === 2)
    assert(alleles.head === ADAMGenotypeAllele.Ref)
    assert(alleles.last === ADAMGenotypeAllele.Alt)
  }
}
