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
import org.bdgenomics.adam.models.{ ReferencePosition, VariantContext }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.{ Contig, Genotype, GenotypeAllele, Variant }
import scala.annotation.tailrec
import scala.math.pow

object BiallelicGenotyper extends GenotyperCompanion {

  val genotyperName: String = "BiallelicGenotyper"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Genotyper = {

    new BiallelicGenotyper(config.getInt("ploidy", 2),
      config.getBoolean("useEM", false),
      config.getBoolean("emitGVCF", true))
  }
}

class BiallelicGenotyper(ploidy: Int = 2,
                         useEM: Boolean = false,
                         emitGVCF: Boolean = true) extends Genotyper with Logging {

  val companion: GenotyperCompanion = BiallelicGenotyper

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
                               observations: Iterable[AlleleObservation]): (Iterable[AlleleObservation], Array[Double], Double) = {

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
      observed.map(o => PhredUtils.phredToErrorProbability((o.phred + o.mapq) / 2))
    }

    // compute error observations
    val refBasesEpsilon = epsilon(refBases)
    val alleleBasesEpsilon = epsilon(alleleBases)
    val mismatchBasesEpsilon = epsilon(mismatchBases)

    // calculate genotype likelihoods
    val likelihoods = states.map(g => {
      // contribution of bases that match the reference
      val productMatch = refBasesEpsilon.map(epsilon => (ploidy - g) * epsilon + g * (1 - epsilon))
        .fold(1.0)(_ * _)

      // contribution of bases that do not match the base
      val productMismatch = alleleBasesEpsilon.map(epsilon => (ploidy - g) * (1 - epsilon) + g * epsilon)
        .fold(1.0)(_ * _)

      productMatch * productMismatch / pow(ploidy.toDouble, k.toDouble)
    })

    // calculate the likelihood of an other alternate allele
    val likelihoodOtherAlt = 1.0 - mismatchBasesEpsilon.fold(1.0)(_ * _)

    (observations, likelihoods, likelihoodOtherAlt)
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
               likelihoodOtherAlt: Double,
               singleSample: Boolean): Option[Genotype] = {

    // were we able to make any observations at this site?
    if (observations.size > 0) {
      // find the genotype state with the maximum likelihood
      val maxLikelihoodState = idxOfMax(likelihoods)

      // if we are calling reference and don't want to emit gvcf, and
      // only have a single sample, now is the time to return
      if (!emitGVCF && maxLikelihoodState == 0 && singleSample) {
        None
      } else {
        // generate called state
        val calls = new Array[GenotypeAllele](ploidy)

        (0 until maxLikelihoodState).foreach(i => {
          calls(i) = GenotypeAllele.Ref
        })

        (maxLikelihoodState until ploidy).foreach(i => {
          calls(i) = GenotypeAllele.Alt
        })

        // get alt and ref alleles
        val alt = variant.getAlternateAllele.toString
        val ref = variant.getReferenceAllele.toString

        // build and return genotype record - just simple statistics for now
        Some(Genotype.newBuilder()
          .setVariant(variant)
          .setSampleId(sampleId)
          .setReadDepth(observations.size)
          .setAlleles(calls.toList)
          .setGenotypeLikelihoods(likelihoods.map(PhredUtils.successProbabilityToPhred).toList)
          .setReferenceReadDepth(observations.filter(_.allele == ref).size)
          .setAlternateReadDepth(observations.filter(_.allele == alt).size)
          .build())
      }
    } else {
      // emit no call
      val calls = new Array[GenotypeAllele](ploidy)

      (0 until ploidy).foreach(i => {
        calls(i) = GenotypeAllele.NoCall
      })

      Some(Genotype.newBuilder()
        .setVariant(variant)
        .setSampleId(sampleId)
        .setReadDepth(0)
        .setAlleles(calls.toList)
        .build())
    }
  }

  def genotypeSite(site: (ReferencePosition, Iterable[Observation])): Option[VariantContext] = {
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
    val compensatedLikelihoodsPerSample = if (useEM) {
      // TODO: connect up EM algorithm
      ???
    } else {
      likelihoodsPerSample
    }

    // construct variant
    val variant = Variant.newBuilder()
      .setReferenceAllele(reference)
      .setAlternateAllele(allele)
      .setContig(Contig.newBuilder()
        .setContigName(pos.referenceName)
        .build())
      .setStart(pos.pos)
      .setEnd(pos.pos + reference.length)
      .build()

    // do we only have a single sample?
    val singleSample = likelihoodsPerSample.size == 1

    // emit calls
    val genotypes = compensatedLikelihoodsPerSample.flatMap(t => {
      // extract info
      val (observations, likelihoods, likelihoodOtherAlt) = t

      // emit the genotypes
      emitCall(variant,
        observations.head.sample,
        observations,
        likelihoods,
        likelihoodOtherAlt,
        singleSample)
    })

    // emit the variant context
    Some(VariantContext(variant, genotypes, None))
  }

  def genotype(observations: RDD[Observation]): RDD[VariantContext] = {
    observations.groupBy(_.pos)
      .flatMap(genotypeSite)
  }
}
