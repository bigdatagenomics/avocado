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
package org.bdgenomics.avocado.genotyping

import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.Logging
import org.bdgenomics.adam.models._
import org.bdgenomics.avocado.algorithms.mutect._
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.{ Contig, Genotype, Variant }

object MutectGenotyper extends GenotyperCompanion {

  val genotyperName: String = "MutectGenotyper"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Genotyper = {

    require(config.containsKey("normalId"),
      "Normal sample ID is not defined in configuration file.")
    require(config.containsKey("somaticId"),
      "Somatic sample ID is not defined in configuration file.")

    // pull out algorithm parameters and return genotyper
    new MutectGenotyper(config.getString("normalId"),
      config.getString("somaticId"),
      stats.contigLengths,
      config.getDouble("threshold", 6.3),
      config.getDouble("somDbSnpThreshold", 5.5),
      config.getDouble("somNovelThreshold", 2.2),
      config.getInt("maxGapEvents", 3),
      None)
  }
}

/**
 *
 * @param threshold Default value of 6.3 corresponds to a 3x10-6 mutation rate,
 *                  see the Online Methods of the Mutect paper for more details.
 * @param somDbSnpThreshold if the site is a known dbSnp site, apply this cutoff for somatic classification
 * @param somNovelThreshold if the site is a novel variant, apply this cutoff for somatic classification
 *
 */
class MutectGenotyper(normalId: String,
                      somaticId: String,
                      val contigLengths: Map[String, Long] = Map.empty,
                      threshold: Double = 6.3,
                      somDbSnpThreshold: Double = 5.5,
                      somNovelThreshold: Double = 2.2,
                      maxGapEventsThreshold: Int = 3,
                      f: Option[Double] = None) extends SiteGenotyper with Logging {

  val companion: GenotyperCompanion = MutectGenotyper

  val model = MutectLogOdds
  // TODO: do we want to do somatic classification here?
  val somaticModel = MutectSomaticLogOdds

  // TODO make sure we only consider the tumor AlleleObservations for this calculation
  def constructVariant(region: ReferenceRegion,
                       ref: String,
                       alt: String,
                       obs: Iterable[AlleleObservation]): VariantContext = {

    val contig = Contig.newBuilder()
      .setContigName(region.referenceName)
      .build()

    val variant = Variant.newBuilder()
      .setStart(region.start)
      .setContig(contig)
      .setEnd(region.end)
      .setReferenceAllele(ref)
      .setAlternateAllele(alt)
      .build()

    val genotypes = Seq()

    VariantContext(variant, genotypes, None)
  }

  // TODO make sure we only consider the tumor AlleleObservations for this calculation
  protected[genotyping] def genotypeSite(region: ReferenceRegion,
                                         referenceObservation: Observation,
                                         alleleObservation: Iterable[AlleleObservation]): Option[VariantContext] = {

    val ref = referenceObservation.allele

    // get all possible alleles for this mutation call
    val alleles: Set[String] = Set(alleleObservation.map(_.allele).toSeq: _*)

    if (ref.size == 1) {
      //TODO do we want to split up normal/tumor here, or earlier in code?
      val tumors_raw = alleleObservation.filter(a => a.sample == somaticId)
      val normals = alleleObservation.filter(a => a.sample == normalId)

      // apply 3 stringent filters to the tumor alleles
      val tumors = tumors_raw.filterNot(ao => {
        val clippedFilter = (ao.clippedBasesReadStart + ao.clippedBasesReadEnd) / ao.unclippedReadLen.toDouble >= 0.3
        val noisyFilter = ao.mismatchQScoreSum >= 100

        // TODO add in mate-rescue filter
        val mateRescueFilter = false
        clippedFilter || noisyFilter || mateRescueFilter
      })

      val rankedAlts: Seq[(Double, String)] =
        (alleles - ref).map { alt =>
          (model.logOdds(ref, alt, alleleObservation, f), alt)
        }.toSeq.sorted.reverse

      //TODO here we should check/flag if there are multiple variants that pass the threshold, this is a
      // downstream filtering criteria we have easy access to right here.
      val passingOddsAlts = rankedAlts.filter(oa => oa._1 >= threshold)

      // TODO do we want to output any kind of info for sites where multiple passing variants are found
      if (passingOddsAlts.size == 1) {
        val alt = passingOddsAlts(0)._2
        // TODO classify somatic status here? or as a downstream step?

        val normalNotHet = somaticModel.logOdds(ref, alt, normals, None)
        val dbSNPsite = false //TODO figure out if this is a dbSNP position
        val passSomatic: Boolean = (dbSNPsite && normalNotHet >= somDbSnpThreshold) || (!dbSNPsite && normalNotHet >= somNovelThreshold)
        val nInsertions = tumors.map(ao => if (ao.distanceToNearestReadInsertion.getOrElse(Int.MaxValue) <= 5) 1 else 0).sum
        val nDeletions = tumors.map(ao => if (ao.distanceToNearestReadDeletion.getOrElse(Int.MaxValue) <= 5) 1 else 0).sum

        val passIndel: Boolean = nInsertions < maxGapEventsThreshold && nDeletions < maxGapEventsThreshold

        val passStringentFilters = tumors.size.toDouble / tumors_raw.size.toDouble > (1.0 - 0.3)

        val passMapq0Filter = tumors_raw.filter(_.mapq.getOrElse(0) == 0).size.toDouble / tumors_raw.size.toDouble <= 0.5 &&
          normals.filter(_.mapq.getOrElse(0) == 0).size.toDouble / normals.size.toDouble <= 0.5

        val onlyTumorMut = tumors.filter(_.allele == alt)

        val passMaxMapqAlt = if (onlyTumorMut.size > 0) onlyTumorMut.map(_.phred).max >= 20 else false

        val passMaxNormalSupport = normals.filter(_.allele == alt).size.toDouble / normals.size.toDouble <= 0.015 ||
          normals.filter(_.allele == alt).map(_.phred).sum < 20

        // TODO implement power-based strand-bias detection
        /* The power to detect a mutant is a function of depth, and the mutant allele fraction (unstranded).
        Basically you assume that the probability of observing a base error is uniform and 0.001 (phred score of 30).
        You see how many reads you require to pass the LOD threshold of 2.0, and then you calculate the binomial
        probability of observing that many or more reads would be observed given the allele fraction and depth.
        You also correct for the fact that the real result is somewhere between the minimum integer number to pass,
        and the number below it, so you scale your probability at k by 1 - (2.0 - lod_(k-1) )/(lod_(k) - lod_(k-1)).
         */
        // TODO implement read-end mutation position clustering detection
        /*
          Median Absolute Deviation, and the raw median of the mutant allele location within reads cannot cluster
          too close to either the beginning or end of each read. Basically you calculate the median of the allele
          supporting reads, and the MAD of those. If the MAD is <= 3 and the median is <= 10, dump it. Similarly
          recalculate these numbers, but using the distance of each allele from the last base of each read, and
          apply the same filters. If either the forward-strand method, or reverse strand method gives you bad
          results, it fails.
         */

        // Do all filters pass?
        if (passSomatic && passIndel && passStringentFilters && passMapq0Filter && passMaxMapqAlt && passMaxNormalSupport) {
          Option(constructVariant(region, ref, alt, alleleObservation))
        } else {
          None
        }

      } else {
        // either there are 0 passing variants, or there are > 1 passing variants
        None
      }
    } else {
      log.info("Dropping site %s, as reference allele is an insertion or complex variant.".format(referenceObservation.pos))
      None
    }
  }

}
