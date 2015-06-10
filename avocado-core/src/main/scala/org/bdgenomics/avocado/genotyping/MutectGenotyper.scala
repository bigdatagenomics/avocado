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
import org.bdgenomics.avocado.algorithms.math.{ LogUtils, LogBinomial }

import scala.annotation.tailrec

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
      config.getDouble("minPassStringentFiltersTumor", 0.3),
      config.getDouble("maxMapq0Fraction", 0.5),
      config.getInt("minPhredSupportingMutant", 20),
      config.getInt("indelNearnessThreshold", 5),
      config.getInt("maxPhredSumMismatchingBases", 100),
      config.getDouble("maxFractionBasesSoftClippedTumor", 0.3),
      config.getDouble("maxNormalSupportingFracToTriggerQscoreCheck", 0.015),
      config.getInt("maxNormalQscoreSumSupportingMutant", 20),
      config.getInt("minMedianDistanceFromReadEnd", 10),
      config.getInt("minMedianAbsoluteDeviationOfAlleleInRead", 3),
      config.getBoolean("experimentalMutectIndelDetector", false),
      config.getDouble("errorForPowerCalculations", 0.001),
      config.getInt("minThetaForPowerCalc", 20),
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
                      minPassStringentFiltersTumor: Double = 0.3,
                      maxMapq0Fraction: Double = 0.5,
                      minPhredSupportingMutant: Int = 20,
                      indelNearnessThreshold: Int = 5,
                      maxPhredSumMismatchingBases: Int = 100,
                      maxFractionBasesSoftClippedTumor: Double = 0.3,
                      maxNormalSupportingFracToTriggerQscoreCheck: Double = 0.015,
                      maxNormalQscoreSumSupportingMutant: Int = 20,
                      minMedianDistanceFromReadEnd: Int = 10,
                      minMedianAbsoluteDeviationOfAlleleInRead: Int = 3,
                      experimentalMutectIndelDetector: Boolean = false,
                      errorForPowerCalculations: Double = 0.001,
                      minThetaForPowerCalc: Int = 20,
                      f: Option[Double] = None) extends SiteGenotyper with Logging {

  val companion: GenotyperCompanion = MutectGenotyper

  val model = MutectLogOdds
  val somaticModel = MutectSomaticLogOdds

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

  protected[genotyping] def genotypeSite(region: ReferenceRegion,
                                         referenceObservation: Observation,
                                         alleleObservation: Iterable[AlleleObservation]): Option[VariantContext] = {

    require(alleleObservation.forall(_.mismatchQScoreSum.isDefined),
      "MD tags must be set in reads for MuTect to function properly. See `samtools calmd` for example.")

    val ref = referenceObservation.allele

    // get all possible alleles for this mutation call
    val alleles: Set[String] = if (experimentalMutectIndelDetector)
      Set(alleleObservation.map(_.allele).toSeq: _*)
    else
      Set(alleleObservation.map(_.allele).toSeq: _*).filter(_.length == 1) // only accept length 1 alleles
    val pointMutation: Boolean = ref.size == 1 && alleles.exists(_.length == 1)

    if (experimentalMutectIndelDetector || pointMutation) {

      val tumorsRaw = alleleObservation.filter(a => a.sample == somaticId)
      val normals = alleleObservation.filter(a => a.sample == normalId)

      // apply 3 stringent filters to the tumor alleles
      val tumors = tumorsRaw.filterNot(ao => {
        val clippedFilter = (ao.clippedBasesReadStart + ao.clippedBasesReadEnd) /
          ao.unclippedReadLen.toDouble >= maxFractionBasesSoftClippedTumor
        val noisyFilter = ao.mismatchQScoreSum.get >= maxPhredSumMismatchingBases
        val mateRescueFilter = ao.mateRescue
        clippedFilter || noisyFilter || mateRescueFilter
      })

      val rankedAlts: Seq[(Double, String)] =
        (alleles - ref).map { alt =>
          (model.logOdds(ref, alt, alleleObservation, f), alt)
        }.toSeq.sorted.reverse

      val passingOddsAlts = rankedAlts.filter(oa => oa._1 >= threshold)

      if (passingOddsAlts.size == 1) {
        val alt = passingOddsAlts(0)._2

        val normalNotHet = somaticModel.logOdds(ref, alt, normals, None)

        val dbSNPsite = false //TODO figure out if this is a dbSNP position

        val passSomatic: Boolean = (dbSNPsite && normalNotHet >= somDbSnpThreshold) || (!dbSNPsite && normalNotHet >= somNovelThreshold)
        val nInsertions = tumors.map(ao => if (math.abs(ao.distanceToNearestReadInsertion.getOrElse(Int.MaxValue)) <= indelNearnessThreshold) 1 else 0).sum
        val nDeletions = tumors.map(ao => if (math.abs(ao.distanceToNearestReadDeletion.getOrElse(Int.MaxValue)) <= indelNearnessThreshold) 1 else 0).sum

        val passIndel: Boolean = nInsertions < maxGapEventsThreshold && nDeletions < maxGapEventsThreshold && pointMutation

        val passStringentFilters = tumors.size.toDouble / tumorsRaw.size.toDouble > (1.0 - minPassStringentFiltersTumor)

        val passMapq0Filter = tumorsRaw.filter(_.mapq.getOrElse(0) == 0).size.toDouble / tumorsRaw.size.toDouble <= maxMapq0Fraction &&
          normals.filter(_.mapq.getOrElse(0) == 0).size.toDouble / normals.size.toDouble <= maxMapq0Fraction

        val onlyTumorMut = tumors.filter(_.allele == alt)

        val passMaxMapqAlt = if (onlyTumorMut.size > 0) onlyTumorMut.map(_.phred).max >= minPhredSupportingMutant else false

        val passMaxNormalSupport = normals.filter(_.allele == alt).size.toDouble / normals.size.toDouble <= maxNormalSupportingFracToTriggerQscoreCheck ||
          normals.filter(_.allele == alt).map(_.phred).sum < maxNormalQscoreSumSupportingMutant

        val tumorPos = tumors.filterNot(_.onNegativeStrand)
        val tumorPosAlt = tumorPos.filter(_.allele == alt)
        val tumorNeg = tumors.filter(_.onNegativeStrand)
        val tumorNegAlt = tumorNeg.filter(_.allele == alt)

        val alleleFrac = onlyTumorMut.size.toDouble / tumors.size.toDouble

        val tumorPosDepth = tumorPos.size
        val tumorNegDepth = tumorNeg.size
        val tPosFrac = if (tumorPosDepth > 0) tumorPosAlt.size.toDouble / tumorPosDepth.toDouble else 0.0
        val tNegFrac = if (tumorNegDepth > 0) tumorNegAlt.size.toDouble / tumorNegDepth.toDouble else 0.0

        val lodPos = model.logOdds(ref, alt, tumorPos, Some(tPosFrac))
        val lodNeg = model.logOdds(ref, alt, tumorNeg, Some(tNegFrac))

        val powerPos = calculateStrandPower(tumorPosDepth, alleleFrac)
        val powerNeg = calculateStrandPower(tumorNegDepth, alleleFrac)

        val passingStrandBias = (powerPos < 0.9 || lodPos >= minThetaForPowerCalc) &&
          (powerNeg < 0.9 || lodNeg >= minThetaForPowerCalc)

        // Only pass mutations that do not cluster at the ends of reads
        val passEndClustering = if (onlyTumorMut.size > 0) {
          val forwardPositions: Seq[Double] = onlyTumorMut.map(_.offsetInAlignment.toDouble).toSeq
          val reversePositions: Seq[Double] = onlyTumorMut.map(ao => ao.alignedReadLen - ao.offsetInAlignment.toDouble - 1.0).toSeq

          val forwardMedian = median(forwardPositions)
          val reverseMedian = median(reversePositions)

          val forwardMad = mad(forwardPositions, forwardMedian)
          val reverseMad = mad(reversePositions, reverseMedian)

          (forwardMad > minMedianAbsoluteDeviationOfAlleleInRead || forwardMedian > minMedianAbsoluteDeviationOfAlleleInRead) &&
            (reverseMad > minMedianAbsoluteDeviationOfAlleleInRead || reverseMedian > minMedianAbsoluteDeviationOfAlleleInRead)

        } else false

        // Do all filters pass?
        if (passSomatic && passIndel && passStringentFilters && passMapq0Filter &&
          passMaxMapqAlt && passMaxNormalSupport && passEndClustering && passingStrandBias) {
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

  def calculateStrandPower(depth: Int, f: Double): Double = {
    /* The power to detect a mutant is a function of depth, and the mutant allele fraction (unstranded).
        Basically you assume that the probability of observing a base error is uniform and 0.001 (phred score of 30).
        You see how many reads you require to pass the LOD threshold of 2.0, and then you calculate the binomial
        probability of observing that many or more reads would be observed given the allele fraction and depth.
        You also correct for the fact that the real result is somewhere between the minimum integer number to pass,
        and the number below it, so you scale your probability at k by 1 - (2.0 - lod_(k-1) )/(lod_(k) - lod_(k-1)).
         */
    if (depth < 1 || f <= 0.0000001) {
      0.0
    } else {

      def getLod(k: Int): Double = {
        val nref = depth - k
        val tf = k / depth.toDouble
        val prefk = nref.toDouble * math.log10(tf * errorForPowerCalculations + (1.0 - tf) * (1.0 - errorForPowerCalculations))
        val paltk = k.toDouble * math.log10(tf * (1.0 - errorForPowerCalculations) + (1.0 - tf) * errorForPowerCalculations)
        val pkm = prefk + paltk
        val pref0 = nref.toDouble * math.log10(1.0 - errorForPowerCalculations)
        val palt0 = k.toDouble * math.log10(errorForPowerCalculations)
        val p0 = pref0 + palt0
        pkm - p0
      }

      @tailrec
      def findMinPassingK(begin: Int, end: Int, passingK: Option[(Int, Double)]): Option[(Int, Double)] = {
        if (begin >= end) passingK
        else {
          val mid = (begin + end) / 2
          val lod = getLod(mid)
          val passingKupdate = if (lod >= minThetaForPowerCalc) Some((mid, lod)) else passingK
          if (lod >= minThetaForPowerCalc && begin < end - 1) findMinPassingK(begin, mid, passingKupdate)
          else if (begin < end - 1) findMinPassingK(mid, end, passingKupdate)
          else passingKupdate
        }
      }

      val kLodOpt = findMinPassingK(1, depth + 1, None)

      if (kLodOpt.isDefined) {
        val (k, lod) = kLodOpt.get
        val probabilities: Array[Double] = LogBinomial.calculateLogProbabilities(math.log(f), depth)
        val binomials = probabilities.drop(k)

        val lodM1 = getLod(k - 1)

        val partialLogPkM1: Double = probabilities(k - 1) + math.log(1.0 - (minThetaForPowerCalc - lodM1) / (lod - lodM1))

        math.exp(LogUtils.sumLogProbabilities(partialLogPkM1 +: binomials))

      } else {
        0.0
      }
    }

  }

  def median(s: Seq[Double]): Double = {
    val (lower, upper) = s.sortWith(_ < _).splitAt(s.size / 2)
    if (s.size % 2 == 0) (lower.last + upper.head) / 2.0 else upper.head
  }

  def mad(s: Seq[Double], m: Double): Double = {
    median(s.map(i => math.abs(i - m)))
  }

}
