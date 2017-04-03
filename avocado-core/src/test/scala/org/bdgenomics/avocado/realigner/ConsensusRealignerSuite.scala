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
package org.bdgenomics.avocado.realigner

import org.bdgenomics.adam.models.SequenceDictionary
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variant.VariantRDD
import org.bdgenomics.formats.avro.Variant

class ConsensusRealignerSuite extends SparkRealignerSuite {

  val allowLegacyCigars = true

  def realign(rdd: AlignmentRecordRDD,
              kmerLength: Int): AlignmentRecordRDD = {
    ConsensusRealigner.realign(rdd, kmerLength)
  }

  val snv = Variant.newBuilder
    .setStart(0L)
    .setReferenceAllele("A")
    .setAlternateAllele("G")
    .build
  val mnv = Variant.newBuilder
    .setStart(1L)
    .setReferenceAllele("AC")
    .setAlternateAllele("GT")
    .build
  val ins = Variant.newBuilder
    .setStart(2L)
    .setReferenceAllele("A")
    .setAlternateAllele("ACC")
    .build
  val del = Variant.newBuilder
    .setStart(3L)
    .setReferenceAllele("GCT")
    .setAlternateAllele("G")
    .build

  test("filter out a snv") {
    assert(!ConsensusRealigner.discardSnvs(snv))
  }

  test("filter out a mnv") {
    assert(!ConsensusRealigner.discardSnvs(mnv))
  }

  test("keep an insert") {
    assert(ConsensusRealigner.discardSnvs(ins))
  }

  test("keep a deletion") {
    assert(ConsensusRealigner.discardSnvs(del))
  }

  sparkTest("filter snv/mnv variants out") {
    val vRdd = VariantRDD(sc.parallelize(Seq(snv, mnv, ins, del)),
      SequenceDictionary.empty)
    val filteredVariants = ConsensusRealigner.filterIndels(vRdd)
      .rdd
      .collect

    val variantIds = filteredVariants.map(_.getStart)
      .toSet

    assert(variantIds.size === 2)
    assert(variantIds(2L))
    assert(variantIds(3L))
  }
}
