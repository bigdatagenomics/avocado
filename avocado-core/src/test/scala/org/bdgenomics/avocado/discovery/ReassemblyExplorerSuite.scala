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
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.avocado.genotyping.BiallelicGenotyper
import org.bdgenomics.formats.avro.{
  AlignmentRecord,
  Contig,
  GenotypeType,
  NucleotideContigFragment
}
import scala.collection.immutable.{ NumericRange, SortedMap }
import scala.collection.mutable.ArrayBuffer

class ReassemblyExplorerSuite extends AvocadoFunSuite {

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

  sparkTest("reassemble and call variants on real data") {
    val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058.sam").getFile
    val reads = sc.loadAlignments(path)
    val sd = new AlignmentRecordRDDFunctions(reads).adamGetSequenceDictionary()
    val cl = sd.records.map(r => (r.name, r.length)).toMap
    val bg = new BiallelicGenotyper(sd, cl)

    val reference = getReferenceFromReads(reads.collect.toSeq)

    val re = new ReassemblyExplorer(20,
      sc.parallelize(Seq(reference)),
      sd,
      cl,
      Double.PositiveInfinity,
      -1.0,
      0.0,
      0.0,
      1000,
      200)
    val obs = re.discover(reads)
    val vc = bg.genotype(obs)
      .flatMap(_.genotypes)
      .filter(g => g.getType != GenotypeType.HOM_REF)
      .collect()

    assert(vc.length === 2)
    val locations = vc.map(_.getVariant.getStart)
      .foreach(l => assert(l == 225048 || l == 225057))
    vc.foreach(gt => assert(gt.getType == GenotypeType.HET))
  }

  sparkTest("put reads into graph for larger real dataset") {
    val path = ClassLoader.getSystemClassLoader.getResource("NA12878_snp_A2G_chr20_225058_longer.sam").getFile
    val reads = sc.loadAlignments(path)
    val sd = new AlignmentRecordRDDFunctions(reads).adamGetSequenceDictionary()
    val cl = sd.records.map(r => (r.name, r.length)).toMap
    val bg = new BiallelicGenotyper(sd, cl)

    val reference = getReferenceFromReads(reads.collect.toSeq)

    val re = new ReassemblyExplorer(20,
      sc.parallelize(Seq(reference)),
      sd,
      cl,
      Double.PositiveInfinity,
      -1.0,
      0.0,
      0.0,
      1000,
      200)
    val obs = re.discover(reads)
    val vc = bg.genotype(obs)
      .flatMap(_.genotypes)
      .filter(g => g.getType != GenotypeType.HOM_REF)
      .collect()

    assert(vc.length === 2)
    val locations = vc.map(_.getVariant.getStart)
      .foreach(l => assert(l == 225048 || l == 225057))
    vc.foreach(gt => assert(gt.getType == GenotypeType.HET))
  }

  test("calculate mismatch and clip rate") {
    val ref = "ACACACAC"
    val reads = Iterable(AlignmentRecord.newBuilder()
      .setSequence("ACAGACACTT")
      .setCigar("8M2S")
      .setStart(0L)
      .build())

    val (mismatchRate, clipRate, coverage) = ReassemblyExplorer.calculateMismatchAndClipRate(reads, ref, 0L)

    assert(mismatchRate > 0.099999 && mismatchRate < 0.100001)
    assert(clipRate > 0.199999 && clipRate < 0.200001)
    assert(coverage > 1.249999 && coverage < 1.250001)
  }
}

