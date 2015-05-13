/*
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.avocado.algorithms.mutect

import org.bdgenomics.adam.models.ReferencePosition
import org.bdgenomics.avocado.models.AlleleObservation
import org.bdgenomics.formats.avro.{ Contig, Genotype, Variant }

/**
 *
 * @param threshold Default value of 6.3 corresponds to a 3x10-6 mutation rate,
 *                  see the Online Methods of the Mutect paper for more details.
 */
class Mutect(val threshold: Double = 6.3, val f: Option[Double] = None) {

  val model = MutectLogOdds
  val alleles: Set[String] = Set("A", "T", "G", "C")

  def constructVariant(pos: ReferencePosition,
                       ref: String,
                       alt: String,
                       obs: Iterable[AlleleObservation]): (Variant, Seq[Genotype]) = {

    val contig = Contig.newBuilder()
      .setContigName(pos.referenceName)
      .build()

    val variant = Variant.newBuilder()
      .setStart(pos.pos)
      .setContig(contig)
      .setEnd(pos.pos + ref.length)
      .setReferenceAllele(ref)
      .setAlternateAllele(alt)
      .build()

    val genotypes = Seq()

    (variant, genotypes)
  }

  def detect(pos: ReferencePosition,
             ref: String,
             obs: Seq[AlleleObservation]): Option[(Variant, Seq[Genotype])] = {

    val rankedAlts: Seq[(Double, String)] =
      (alleles - ref).map { alt =>
        (model.logOdds(ref, alt, obs, f), alt)
      }.toSeq.sorted.reverse

    val (logOdds, alt) = rankedAlts(0)

    if (logOdds >= threshold) {
      Option(constructVariant(pos, ref, alt, obs))
    } else {
      None
    }
  }

}
