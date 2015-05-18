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
 *
 */
class Mutect(val threshold: Double = 6.3, val f: Option[Double] = None) {

  val model = MutectLogOdds
  // TODO: do we want to do somatic classification here?
  // val somaticClassifier = MutectSomaticLogOdds
  val alleles: Set[String] = Set("A", "T", "G", "C")

  // TODO make sure we only consider the tumor AlleleObservations for this calculation
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

  // TODO make sure we only consider the tumor AlleleObservations for this calculation
  def detect(pos: ReferencePosition,
             ref: String,
             obs: Seq[AlleleObservation]): Option[(Variant, Seq[Genotype])] = {

    //TODO do we want to split up normal/tumor here, or earlier in code?
    val tumors = obs.filter(a => a.sample == "tumor")
    //val normals = obs.filter(a => a.sample == "normal")

    val rankedAlts: Seq[(Double, String)] =
      (alleles - ref).map { alt =>
        (model.logOdds(ref, alt, obs, f), alt)
      }.toSeq.sorted.reverse

    //TODO here we should check/flag if there are multiple variants that pass the threshold, this is a
    // downstream filtering criteria we have easy access to right here.
    val passingOddsAlts = rankedAlts.filter(oa => oa._1 >= threshold)

    // TODO do we want to output any kind of info for sites where multiple passing variants are found
    if (passingOddsAlts.size == 1) {
      // TODO classify somatic status here? or as a downstream step?

      val alt = passingOddsAlts(0)._2
      Option(constructVariant(pos, ref, alt, obs))
    } else {
      // either there are 0 passing variants, or there are > 1 passing variants
      None
    }
  }

}
