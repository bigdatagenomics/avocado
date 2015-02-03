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

import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext._
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDDFunctions
import org.bdgenomics.adam.rich.RichGenotype._
import org.bdgenomics.adam.rich.{ RichAlignmentRecord, RichGenotype }
import org.bdgenomics.adam.util.SparkFunSuite
import org.bdgenomics.avocado.genotyping.BiallelicGenotyper
import org.bdgenomics.formats.avro.{
  AlignmentRecord,
  Contig,
  GenotypeType,
  NucleotideContigFragment
}
import scala.collection.immutable.{ NumericRange, SortedMap }
import scala.collection.mutable.ArrayBuffer

class ReassemblyExplorerSuite extends SparkFunSuite {

  override val properties = Map(("spark.serializer", "org.apache.spark.serializer.KryoSerializer"),
    ("spark.kryo.registrator", "org.bdgenomics.adam.serialization.ADAMKryoRegistrator"),
    ("spark.kryoserializer.buffer.mb", "128"),
    ("spark.kryo.referenceTracking", "true"))

  // shamelessly borrowed from the indel realigner while we are refactoring...
  def getReferenceFromReads(poorReads: Seq[AlignmentRecord]): NucleotideContigFragment = {
    val reads = poorReads.map(RichAlignmentRecord(_))

    // get reference and range from a single read
    val readRefs = reads.map((r: RichAlignmentRecord) => {
      (r.mdTag.get.getReference(r), r.getStart.toLong to r.getEnd)
    })
      .toSeq
      .sortBy(_._2.head)

    // fold over sequences and append - sequence is sorted at start
    val ref = readRefs.reverse.foldRight[(String, Long)](("", readRefs.head._2.head))((refReads: (String, NumericRange[Long]), reference: (String, Long)) => {
      if (refReads._2.end < reference._2) {
        reference
      } else if (reference._2 >= refReads._2.head) {
        (reference._1 + refReads._1.substring((reference._2 - refReads._2.head).toInt), refReads._2.end)
      } else {
        // there is a gap in the sequence
        throw new IllegalArgumentException("Current sequence has a gap at " + reference._2 + "with " + refReads._2.head + "," + refReads._2.end + ".")
      }
    })

    // get reference info
    val ctg = poorReads.head.getContig
    val start = poorReads.map(_.getStart).min

    NucleotideContigFragment.newBuilder()
      .setContig(ctg)
      .setFragmentSequence(ref._1)
      .setFragmentStartPosition(start)
      .build()
  }

  def na12878_chr20_snp_reads: RDD[AlignmentRecord] = {
    val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058.sam").getFile
    sc.loadAlignments(path)
  }

  def more_na12878_chr20_snp_reads: RDD[AlignmentRecord] = {
    val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058_longer.sam").getFile
    sc.loadAlignments(path)
  }

  lazy val sd = new AlignmentRecordRDDFunctions(na12878_chr20_snp_reads).adamGetSequenceDictionary()
  lazy val cl = sd.records.map(r => (r.name, r.length)).toMap

  sparkTest("reassemble and call variants on real data") {
    val reads = na12878_chr20_snp_reads

    val reference = getReferenceFromReads(reads.collect.toSeq)

    val re = new ReassemblyExplorer(20, sc.parallelize(Seq(reference)), sd, cl)
    val obs = re.discover(reads)
    val bg = new BiallelicGenotyper(sd, 2, false)
    val vc = bg.genotype(obs)
      .flatMap(_.genotypes)
      .filter(g => g.getType != GenotypeType.HOM_REF)

    assert(vc.count === 1)
    val gt = vc.first
    assert(gt.getVariant.getReferenceAllele === "A")
    assert(gt.getVariant.getAlternateAllele === "G")
    assert(gt.getVariant.getStart === 225057L)
    assert(gt.getType == GenotypeType.HET)
  }

  sparkTest("put reads into graph for larger real dataset") {
    val reads = more_na12878_chr20_snp_reads

    val reference = getReferenceFromReads(reads.collect.toSeq)

    val re = new ReassemblyExplorer(20, sc.parallelize(Seq(reference)), sd, cl)
    val obs = re.discover(reads)
    val bg = new BiallelicGenotyper(sd, 2, false)
    val vc = bg.genotype(obs)
      .flatMap(_.genotypes)
      .filter(g => g.getType != GenotypeType.HOM_REF)
      .collect()

    assert(vc.length === 4)
    val locations = vc.map(_.getVariant.getStart)
      .foreach(l => assert(l == 224940 || l == 225048 || l == 224970 || l == 225057))
    vc.foreach(gt => assert(gt.getType == GenotypeType.HET))
  }
}

