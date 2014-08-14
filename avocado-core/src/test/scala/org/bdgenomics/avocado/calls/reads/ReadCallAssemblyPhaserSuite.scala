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
package org.bdgenomics.avocado.calls.reads

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichAlignmentRecord
import org.bdgenomics.adam.util.SparkFunSuite
import org.bdgenomics.formats.avro.{ GenotypeAllele, AlignmentRecord, Contig }
import org.bdgenomics.avocado.algorithms.hmm._
import org.bdgenomics.avocado.partitioners.PartitionSet
import parquet.filter.UnboundRecordFilter
import scala.collection.JavaConversions._
import scala.collection.immutable.SortedMap

class ReadCallAssemblyPhaserSuite extends ReadCallHaplotypesSuite {

  val rc_short = new ReadCallAssemblyPhaser(emptyPartition, 4, 4)
  val rc_long = new ReadCallAssemblyPhaser(emptyPartition, 20, 40)
  val rc_long_trim = new ReadCallAssemblyPhaser(emptyPartition, 20, 40, 5, true, Some(0.05))

  sparkTest("call doesn't change after \"aggressive\" trimming") {
    val reads = na12878_chr20_snp_reads.collect.toSeq
    val reference = rc_long_trim.getReference(reads)

    val kmerGraph = rc_long_trim.generateHaplotypes(reads, reference)

    val variants = rc_long_trim.scoreHaplotypes(reads, kmerGraph, reference)

    assert(variants.length === 1)
    assert(variants.head.position.pos === 225057L)
    assert(variants.head.variant.variant.getReferenceAllele === "A")
    assert(variants.head.variant.variant.getAlternateAllele === "G")
    val alleles: List[GenotypeAllele] = asScalaBuffer(variants.head.genotypes.head.getAlleles).toList
    assert(alleles.length === 2)
    assert(alleles.head === GenotypeAllele.Ref)
    assert(alleles.last === GenotypeAllele.Alt)
  }
}
