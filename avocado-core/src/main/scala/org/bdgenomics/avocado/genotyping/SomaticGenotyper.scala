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
package org.bdgenomics.avocado.genotyping

import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{
  ReferencePosition,
  SequenceDictionary,
  SequenceRecord,
  VariantContext
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.{ Contig, Genotype, GenotypeAllele, Variant }
import scala.annotation.tailrec
import scala.math.pow

object SomaticGenotyper extends GenotyperCompanion {

  val genotyperName: String = "SomaticGenotyper"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Genotyper = {

    new SomaticGenotyper(stats.sequenceDict,
      config.getString("normalSample"),
      config.getString("somaticSample"),
      config.getInt("normalPloidy", 2))
  }
}

class SomaticGenotyper(sd: SequenceDictionary,
                       normalSample: String,
                       somaticSample: String,
                       normalPloidy: Int = 2) extends Genotyper with Logging {

  val companion: GenotyperCompanion = SomaticGenotyper

  protected val bg = new BiallelicGenotyper(sd, normalPloidy)

  def genotypeSite(site: (ReferencePosition, Iterable[Observation])): Iterable[VariantContext] = {
    val (pos, observations) = site

    // get reference allele
    val reference = observations.find(obs => obs match {
      case ao: AlleleObservation => false
      case _                     => true
    }).fold({
      log.info("Had alleles observed, but no reference at " + pos)
      return None
    })(_.allele)

    // get allele observations
    val alleleObservations: Iterable[AlleleObservation] = observations.flatMap(obs => obs match {
      case ao: AlleleObservation => Some(ao)
      case _                     => None
    })

    // split out normal and somatic sample observations
    val normalObservations = alleleObservations.filter(_.sample == normalSample)
    val somaticObservations = alleleObservations.filter(_.sample == somaticSample)

    // get most frequently observed non-ref allele from normal sample
    val nonRefAlleles = normalObservations.filter(_.allele != reference)
    val allele = if (nonRefAlleles.size > 0) {
      nonRefAlleles.groupBy(_.allele)
        .maxBy(kv => kv._2.size)
        ._1
    } else {
      "N" // we need a non-reference allele for calculating genotype likelihoods
    }

    // calculate site genotype likelihoods for the normal sample
    val (_, normalLikelihoods, _) = bg.scoreGenotypeLikelihoods(reference,
      allele,
      normalObservations)
    val normalGenotype = bg.idxOfMax(normalLikelihoods)

    // calculate somatic depth
    val somaticDepth = somaticObservations.size

    // find alleles in somatic sample
    val somaticAlleleObservations = somaticObservations.filter(obs => {
      !(normalGenotype != 0 && obs.allele == reference) &&
        !(normalGenotype != 2 && obs.allele == allele)
    }).groupBy(_.allele)

    // emit calls
    somaticAlleleObservations.map(kv => {
      val (allele, obs) = kv

      // emit variant
      val variant = Variant.newBuilder()
        .setReferenceAllele(reference)
        .setAlternateAllele(allele)
        .setContig(SequenceRecord.toADAMContig(sd(pos.referenceName).get))
        .setStart(pos.pos)
        .setEnd(pos.pos + reference.length)
        .build()

      // emit genotype - not sure what to do about genotype state...
      val genotype = Genotype.newBuilder()
        .setVariant(variant)
        .setSampleId(somaticSample)
        .setReadDepth(somaticDepth)
        .setAlternateReadDepth(obs.size)
        .setGenotypeQuality(obs.map(o => (o.phred + o.mapq) / 2).sum)
        .build()

      VariantContext(variant, Iterable(genotype), None)
    })
  }

  def genotype(observations: RDD[Observation]): RDD[VariantContext] = {
    observations.groupBy(_.pos)
      .flatMap(genotypeSite)
  }
}
