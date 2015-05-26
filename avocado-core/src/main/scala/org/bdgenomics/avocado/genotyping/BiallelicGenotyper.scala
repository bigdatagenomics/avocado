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
import org.apache.spark.SparkContext._
import org.apache.spark.Logging
import org.apache.spark.rdd.{ InstrumentedRDD, RDD }
import org.bdgenomics.adam.models.{
  ReferencePosition,
  ReferenceRegion,
  SequenceDictionary,
  SequenceRecord,
  VariantContext
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.GenomicPositionPartitioner
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.algorithms.em.EMForAlleles
import org.bdgenomics.avocado.algorithms.math._
import org.bdgenomics.avocado.genotyping.annotators.VariantCallingAnnotator
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.avocado.algorithms.math.LogToPhred.log2phred
import org.bdgenomics.formats.avro.{ Contig, Genotype, GenotypeAllele, Variant, VariantCallingAnnotations }
import scala.annotation.tailrec
import scala.math.{ exp, expm1, log => mathLog, log1p, max, min, pow }

object BiallelicGenotyper extends GenotyperCompanion {

  val genotyperName: String = "BiallelicGenotyper"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Genotyper = {

    // get finishing conditions for EM algorithms
    val useEM = config.getBoolean("useEM", true)
    val maxIterations = if (config.containsKey("maxEMIterations")) {
      val iterLimit = config.getInt("maxEMIterations")
      if (iterLimit <= 0) {
        throw new IllegalArgumentException("EM iteration limit must be greater than 0.")
      }
      Some(iterLimit)
    } else {
      None
    }
    val tolerance = if (config.containsKey("emTolerance")) {
      val emTol = config.getDouble("emTolerance")
      if (emTol < 0.0 || emTol > 1.0) {
        throw new IllegalArgumentException("EM tolerance must be between 0 and 1, non-inclusive.")
      }
      Some(emTol)
    } else {
      None
    }

    if (maxIterations.isEmpty && tolerance.isEmpty && useEM) {
      throw new IllegalArgumentException("At least one constraint must be defined for the EM algorithm.")
    }

    // what level do we saturate the reference frequency to if we encounter underflow?
    val referenceFrequency = config.getDouble("referenceFrequency", 0.999)
    if (referenceFrequency < 0.0 || referenceFrequency > 1.0) {
      throw new IllegalArgumentException("Reference frequency must be between 0 and 1, non-inclusive.")
    }
    val saturationThreshold = config.getDouble("emSaturationThreshold", 0.001)
    if (saturationThreshold < 0.0 || saturationThreshold > 1.0) {
      throw new IllegalArgumentException("Saturation threshold must be between 0 and 1, non-inclusive.")
    }

    new BiallelicGenotyper(stats.sequenceDict,
      stats.contigLengths,
      config.getInt("ploidy", 2),
      useEM,
      config.getBoolean("emitGVCF", true),
      referenceFrequency,
      maxIterations,
      tolerance,
      saturationThreshold)
  }
}

class BiallelicGenotyper(sd: SequenceDictionary,
                         contigLengths: Map[String, Long],
                         ploidy: Int = 2,
                         useEM: Boolean = false,
                         emitGVCF: Boolean = true,
                         estimatedReferenceFrequency: Double = 0.999,
                         maxIterations: Option[Int] = Some(10),
                         tolerance: Option[Double] = Some(1e-3),
                         saturationThreshold: Double = 0.001) extends Genotyper with Logging {

  val companion: GenotyperCompanion = BiallelicGenotyper

  // necessary for conversion to/from phred

  // precompute log ploidy
  private val logPloidy = mathLog(ploidy.toDouble)

  /**
   * Scores likelihoods for genotypes using equation from Li 2009.
   *
   * @param reference Reference base at this site.
   * @param allele Alternate allele at this site.
   * @param observations Observed alleles at this site in this sample.
   * @return List of doubles corresponding to likelihood of homozygous ref (0),
   *         heterozygous (1), and homozygous non-reference (2).
   */
  def scoreGenotypeLikelihoods(reference: String,
                               allele: String,
                               observations: Iterable[AlleleObservation]): (Iterable[AlleleObservation], Array[Double], Array[Double]) = ScoringLikelihoods.time {

    // count bases observed
    val k = observations.size

    // genotype states
    val states = (0 to ploidy).toArray

    /* find genotype for pileup
     * likelihood for genotype is derived from:
     * L(g) = 1/m^k *
     *        product over i->1..l         (g * e_i       + (m - g) * (1 - e_i)) *
     *        product over i->l+1..(k - j) (g * (1 - e_i) + (m - g) * e_i)
     *        // ignore for moment: product over i->(k - j)..k (e_i)
     *
     * Where:
     * - m is ploidy
     * - g in 0..m is genotype -> g ==  # of reference bases in allele
     * - e_i is error probablity of nucleotide i
     * - k is the number of observations
     * - l is the number of bases that match the reference
     * - j is the number of bases that mismatch the reference and the allele
     */
    val refBases = observations.filter(o => o.allele == reference)
    val alleleBases = observations.filter(o => o.allele == allele)
    val mismatchBases = observations.filter(o => o.allele != reference && o.allele != allele)

    def epsilon(observed: Iterable[AlleleObservation]): Iterable[Double] = {
      // we must take the min here; if we don't, phred=0 mapq=0 bases have epsilon of 1.0
      // which leads to log(1.0 - 1.0) = log(0.0) = -Infinity, which messes everything up
      observed.map(o => min(0.99999, o.mapq.fold(PhredUtils.phredToErrorProbability(o.phred))(mq =>
        PhredUtils.phredToErrorProbability((o.phred + mq) >> 1))))
    }

    // compute error observations
    val refBasesEpsilon = epsilon(refBases)
    val alleleBasesEpsilon = epsilon(alleleBases)
    val mismatchBasesEpsilon = epsilon(mismatchBases)
    val nonRefBasesEpsilon = alleleBasesEpsilon ++ mismatchBasesEpsilon

    def calculateLikelihoods(refEpsilon: Iterable[Double],
                             alleleEpsilon: Iterable[Double]): Array[Double] = {
      states.map(g => {
        // contribution of bases that match the reference
        val productMatch = refEpsilon.map(epsilon => mathLog((ploidy - g) * epsilon + g * (1 - epsilon)))
          .sum

        // contribution of bases that do not match the base
        val productMismatch = alleleEpsilon.map(epsilon => mathLog((ploidy - g) * (1 - epsilon) + g * epsilon))
          .sum

        productMatch + productMismatch - (logPloidy * k.toDouble)
      })
    }

    // calculate genotype likelihoods
    val likelihoods = calculateLikelihoods(refBasesEpsilon, alleleBasesEpsilon)
    val anyAltLikelihoods = calculateLikelihoods(refBasesEpsilon, nonRefBasesEpsilon)

    (observations, likelihoods, anyAltLikelihoods)
  }

  @tailrec final def idxOfMax(array: Array[Double],
                              idx: Int = 1,
                              maxIdx: Int = 0): Int = {
    // are we at the end of the array? if so, return.
    if (idx >= array.length) {
      maxIdx
    } else {
      // do we have a new max? if so, update the current max index.
      val newMaxIdx = if (array(idx) > array(maxIdx)) {
        idx
      } else {
        maxIdx
      }

      // recurse
      idxOfMax(array, idx + 1, newMaxIdx)
    }
  }

  def emitCall(variant: Variant,
               sampleId: String,
               observations: Iterable[AlleleObservation],
               likelihoods: Array[Double],
               anyAltLikelihoods: Array[Double],
               priors: Array[Double]): Genotype = EmittingCall.time {

    // were we able to make any observations at this site?
    if (observations.size > 0) {
      // compute and marginalize posterior
      val posteriors = new Array[Double](likelihoods.length)
      (0 until posteriors.length).foreach(i => {
        posteriors(i) = likelihoods(i) + priors(i)
      })
      val norm = LogUtils.sumLogProbabilities(posteriors)
      (0 until posteriors.length).foreach(i => {
        posteriors(i) = posteriors(i) - norm
      })

      // find the genotype state with the maximum posterior probability
      val maxPosteriorState = idxOfMax(posteriors)
      val genotypeProbability = posteriors(maxPosteriorState)

      // generate called state
      val calls = new Array[GenotypeAllele](ploidy)

      (0 until maxPosteriorState).foreach(i => {
        calls(i) = GenotypeAllele.Ref
      })

      (maxPosteriorState until ploidy).foreach(i => {
        calls(i) = GenotypeAllele.Alt
      })

      // get alt and ref alleles
      val alt = variant.getAlternateAllele.toString
      val ref = variant.getReferenceAllele.toString

      // get variant annotations
      val variantAnnotations = VariantCallingAnnotator(variant,
        observations,
        log2phred(genotypeProbability),
        likelihoods).build()

      // build and return genotype record - just simple statistics for now
      Genotype.newBuilder()
        .setVariant(variant)
        .setSampleId(sampleId)
        .setReadDepth(observations.size)
        .setAlleles(calls.toList)
        .setGenotypeLikelihoods(likelihoods.map(v => log2phred(v).toInt).toList)
        .setNonReferenceLikelihoods(anyAltLikelihoods.map(v => log2phred(v).toInt).toList)
        .setGenotypeQuality(log2phred(genotypeProbability).toInt)
        .setReferenceReadDepth(observations.filter(_.allele == ref).size)
        .setAlternateReadDepth(observations.filter(_.allele == alt).size)
        .setVariantCallingAnnotations(variantAnnotations)
        .build()

    } else {
      // emit no call
      val calls = new Array[GenotypeAllele](ploidy)

      (0 until ploidy).foreach(i => {
        calls(i) = GenotypeAllele.NoCall
      })

      Genotype.newBuilder()
        .setVariant(variant)
        .setSampleId(sampleId)
        .setReadDepth(0)
        .setAlleles(calls.toList)
        .build()
    }
  }

  def genotypeSite(region: ReferenceRegion,
                   referenceObservation: Observation,
                   alleleObservations: Iterable[AlleleObservation]): VariantContext = GenotypingSite.time {

    // get reference allele
    val reference = referenceObservation.allele

    // get observations per sample
    val observationsBySample = alleleObservations.groupBy(_.sample)

    // find most frequently observed non-ref allele
    val nonRefAlleles = alleleObservations.filter(_.allele != reference)

    val allele = if (nonRefAlleles.size > 0) {
      nonRefAlleles.groupBy(_.allele)
        .maxBy(kv => kv._2.size)
        ._1
    } else {
      "N" // we need a non-reference allele for calculating genotype likelihoods
    }

    // generate genotype likelihoods
    val likelihoodsPerSample = observationsBySample.map(kv => {
      scoreGenotypeLikelihoods(reference, allele, kv._2)
    })

    // compensate likelihoods on the basis of population statistics
    val logMajorAlleleFrequency = if (useEM) {
      min(mathLog(1.0 - saturationThreshold),
        max(EMForAlleles.emForMAF(likelihoodsPerSample.flatMap(s => {
          // did we have any observations from this sample?
          if (s._1.size > 0) {
            Some(s._2)
          } else {
            None
          }
        }).toArray,
          estimatedReferenceFrequency,
          maxIterations,
          tolerance), mathLog(saturationThreshold)))
    } else {
      mathLog(estimatedReferenceFrequency)
    }
    val statePriors = LogBinomial.calculateLogProbabilities(logMajorAlleleFrequency, ploidy)

    // construct variant
    val variant = Variant.newBuilder()
      .setReferenceAllele(reference)
      .setAlternateAllele(allele)
      .setContig(SequenceRecord.toADAMContig(sd(region.referenceName).get))
      .setStart(region.start)
      .setEnd(region.end)
      .build()

    // emit calls
    val genotypes = likelihoodsPerSample.map(t => {
      // extract info
      val (observations, likelihoods, anyAltLikelihoods) = t

      // emit the genotypes
      emitCall(variant,
        observations.head.sample,
        observations,
        likelihoods,
        anyAltLikelihoods,
        statePriors)
    })

    // emit the variant context
    VariantContext(variant, genotypes, None)
  }

  def genotype(observations: RDD[Observation]): RDD[VariantContext] = {
    val sortedRdd = observations.keyBy(_.pos)
      .repartitionAndSortWithinPartitions(GenomicPositionPartitioner(observations.partitions.size,
        contigLengths))
    val instrumentedRdd = {
      import org.apache.spark.rdd.MetricsContext._
      sortedRdd.instrument
    }
    instrumentedRdd.mapPartitions(genotypeIterator)
  }

  def updateObservations(region: ReferenceRegion,
                         obs: List[Observation]): (Observation, Iterable[AlleleObservation]) = {
    // extract references
    val ref = obs.flatMap(o => o match {
      case ao: AlleleObservation => None
      case oo: Observation       => Some(oo)
    }).sortBy(_.pos)
      .map(_.allele)
      .mkString
    val refLen = ref.length
    val refPos = ReferencePosition(region.referenceName, region.start)

    // put read observations together
    val readObservations = obs.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case oo: Observation       => None
    })

    // if the reference length is greater than 1, we need to expand the values
    val finalObservations = if (refLen > 1) {
      readObservations.groupBy(_.readId)
        .map(ro => {
          // sort the read observations
          val sortedRo = ro._2.sortBy(_.pos)

          // get the start and end position
          val start = sortedRo.head.pos.pos
          val end = sortedRo.last.pos.pos + sortedRo.last.length

          // reduce down to get the sequence
          val sequence = (ref.take((start - region.start).toInt) +
            sortedRo.map(_.allele).mkString +
            ref.takeRight((region.end - end).toInt))
          val phred = sortedRo.map(_.phred).sum / sortedRo.length

          // make observation
          AlleleObservation(refPos,
            refLen,
            sequence,
            phred,
            sortedRo.head.mapq,
            sortedRo.head.onNegativeStrand,
            sortedRo.head.firstOfPair,
            sortedRo.head.offsetInRead,
            sortedRo.head.sample,
            sortedRo.head.readId)
        })
    } else {
      readObservations
    }

    // emit reference observation
    (new Observation(refPos, ref), finalObservations)
  }

  def genotypeIterator(iter: Iterator[(ReferencePosition, Observation)]): Iterator[VariantContext] = {
    // do we have any elements in this iterator?
    if (iter.hasNext) {
      // peek at the head element
      val (headPos, headObs) = iter.next

      // and initialize state
      var windowStart = headPos.pos
      var windowEnd = windowStart + headObs.length
      var obsList = List(headObs)

      def flush(): VariantContext = {
        val region = ReferenceRegion(headPos.referenceName, windowStart, windowEnd)

        // update observations
        val (refObs, obsIterable) = updateObservations(region, obsList)

        // compute genotypes
        genotypeSite(region,
          refObs,
          obsIterable)
      }

      // start processing observation windows!
      iter.flatMap(kv => {
        val (pos, obs) = kv

        // does this observation start past the end of this window? if so, flush our observations
        // and rebuild the window
        if (pos.pos >= windowEnd) {
          // flush the window
          val vc = flush()

          // reinitialize state
          windowStart = pos.pos
          windowEnd = windowStart + obs.length
          obsList = List(obs)

          // emit variant context
          Some(vc)
        } else {
          // make possible update to window end
          windowEnd = max(windowEnd, pos.pos + obs.length)

          // append observation to list
          obsList = obs :: obsList

          // nothing to emit
          None
        }
      }) ++ Iterator(flush())
    } else {
      Iterator()
    }
  }
}
