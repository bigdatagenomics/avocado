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

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.MetricsContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.GenomeBins
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variant.{
  GenotypeRDD,
  VariantRDD
}
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.models.Observation
import org.bdgenomics.avocado.util.{
  Downsampler,
  HardLimiter,
  LogPhred,
  LogUtils
}
import org.bdgenomics.formats.avro.{
  AlignmentRecord,
  Genotype,
  GenotypeAllele,
  Sample,
  Variant,
  VariantCallingAnnotations
}
import org.bdgenomics.utils.instrumentation.Metrics
import org.bdgenomics.utils.misc.{ Logging, MathUtils }
import scala.annotation.tailrec
import scala.collection.JavaConversions._
import scala.math.{ exp, log => mathLog, log10, sqrt }

/**
 * Calls genotypes from reads assuming a biallelic genotyping model.
 *
 * Uses the genotyping model from:
 *
 *   Li, Heng. "A statistical framework for SNP calling, mutation discovery,
 *   association mapping and population genetical parameter estimation from
 *   sequencing data." Bioinformatics 27.21 (2011): 2987-2993.
 *
 * Assumes no prior in favor of/against the reference, and that there is no
 * copy number variation in the sample (all sites have consistent ploidy).
 */
private[avocado] object BiallelicGenotyper extends Serializable with Logging {

  /**
   * Force calls variants from the input read dataset.
   *
   * @param reads The reads to use as evidence.
   * @param variants The variants to force call.
   * @param ploidy The assumed copy number at each site.
   * @param optDesiredPartitionCount Optional parameter for setting the number
   *   of partitions desired after the shuffle region join. Defaults to None.
   * @param optDesiredPartitionSize An optional number of reads per partition to
   *   target.
   * @param optDesiredMaxCoverage An optional cap for coverage per site.
   * @return Returns genotype calls.
   */
  def call(reads: AlignmentRecordRDD,
           variants: VariantRDD,
           ploidy: Int,
           optDesiredPartitionCount: Option[Int] = None,
           optDesiredPartitionSize: Option[Int] = None,
           optDesiredMaxCoverage: Option[Int] = None): GenotypeRDD = CallGenotypes.time {

    // validate metadata
    require(variants.sequences.isCompatibleWith(reads.sequences),
      "Variant sequence dictionary (%s) is not compatible with read dictionary (%s).".format(
        variants.sequences, reads.sequences))
    val samples = reads.recordGroups.recordGroups.map(_.sample).toSet
    require(samples.size == 1,
      "Currently, we only support a single sample. Saw: %s.".format(
        samples.mkString(", ")))

    // join reads against variants
    val useBroadcastJoin = false
    val useTreeJoin = true
    val joinedRdd = JoinReadsAndVariants.time {
      if (useTreeJoin) {

        variants.broadcastRegionJoinAndGroupByRight(reads)
          .rdd.map(_.swap)
      } else if (useBroadcastJoin) {

        val joinedGRdd = variants.broadcastRegionJoin(reads)
        joinedGRdd.rdd.map(kv => {
          val (variant, read) = kv
          (read, Iterable(variant))
        })
      } else {
        val refLength = reads.sequences.records.map(_.length).sum
        val partitionSize = refLength / optDesiredPartitionCount.getOrElse(reads.rdd.partitions.size)

        val shuffledRdd = reads.shuffleRegionJoinAndGroupByLeft(variants,
          optPartitions = Some(partitionSize.toInt))
          .rdd

        val rdd = if (Metrics.isRecording) shuffledRdd.instrument() else shuffledRdd
        val optCoverageThresholdedRdd = optDesiredMaxCoverage.fold(rdd)(maxCoverage => {
          HardLimiter(rdd, maxCoverage)
        })
        optDesiredPartitionSize.fold(optCoverageThresholdedRdd)(size => {

          val contigLengths = reads.sequences
            .records
            .map(r => (r.name -> r.length))
            .toMap
          val refLength = contigLengths.values.sum
          val partitionSize = refLength / optDesiredPartitionCount.getOrElse(reads.rdd.partitions.size)

          val bins = GenomeBins(partitionSize, contigLengths)

          Downsampler.downsample(optCoverageThresholdedRdd, size, bins)
        })
      }
    }

    // score variants and get observations
    val observationRdd = readsToObservations(joinedRdd, ploidy)

    // genotype observations
    val genotypeRdd = observationsToGenotypes(observationRdd, samples.head)

    GenotypeRDD(genotypeRdd,
      variants.sequences,
      samples.map(s => {
        Sample.newBuilder()
          .setSampleId(s)
          .setName(s)
          .build()
      }).toSeq)
  }

  /**
   * Discovers variants and calls genotypes from the input read dataset.
   *
   * @param reads The reads to use as evidence.
   * @param ploidy The assumed copy number at each site.
   * @param optDesiredPartitionCount Optional parameter for setting the number
   *   of partitions desired after the shuffle region join. Defaults to None.
   * @param optPhredThreshold An optional threshold that discards all variants
   *   not supported by bases of at least a given phred score.
   * @param optDesiredPartitionSize An optional number of reads per partition to
   *   target.
   * @param optDesiredMaxCoverage An optional cap for coverage per site.
   * @return Returns genotype calls.
   */
  def discoverAndCall(reads: AlignmentRecordRDD,
                      ploidy: Int,
                      optDesiredPartitionCount: Option[Int] = None,
                      optPhredThreshold: Option[Int] = None,
                      optMinObservations: Option[Int] = None,
                      optDesiredPartitionSize: Option[Int] = None,
                      optDesiredMaxCoverage: Option[Int] = None): GenotypeRDD = {

    // get rdd storage level and warn if not persisted
    val readSl = reads.rdd.getStorageLevel
    if (!readSl.useDisk || !readSl.useMemory) {
      log.warn("Input RDD is not persisted. Performance may be degraded.")
    }

    // discover variants
    val variants = DiscoverVariants(reads,
      optPhredThreshold = optPhredThreshold,
      optMinObservations = optMinObservations)

    // "force" call the variants we have discovered in the input reads
    call(reads, variants, ploidy,
      optDesiredPartitionCount = optDesiredPartitionCount,
      optDesiredPartitionSize = optDesiredPartitionSize)
  }

  /**
   * Scores the putative variants covered by a read against a single read.
   *
   * @param readAndVariants A tuple pairing a read with the variants it covers.
   * @param copyNumber The number of copies of this locus expected to exist at
   *   this site.
   * @return Returns an iterable collection of (Variant, Observation) pairs.
   */
  private[genotyping] def readToObservations(
    readAndVariants: (AlignmentRecord, Iterable[Variant]),
    copyNumber: Int): Iterable[(Variant, Observation)] = ObserveRead.time {

    // unpack tuple
    val (read, variants) = readAndVariants

    if (variants.isEmpty) {
      Iterable.empty
    } else {
      try {

        // observe this read
        val observations = Observer.observeRead(read, copyNumber)
          .map(kv => {
            val ((rr, allele, _), obs) = kv
            (rr, allele, obs)
          })

        // find all variant/observation intersections
        val intersection = IntersectVariants.time {
          variants.map(v => {
            val rr = ReferenceRegion(v)

            val obs = observations.filter(_._1.overlaps(rr))

            (v, obs)
          })
        }

        // process these intersections
        ProcessIntersections.time {
          val obsMap = intersection.flatMap(kv => {
            val (variant, observed) = kv

            // what type of variant is this?
            // - if snp or deletion, look for single observation matching alt
            // - if insertion, look for observation matching insert tail
            //
            // FIXME: if we don't see the variant, take the first thing and invert it

            if (observed.isEmpty) {
              None
            } else if (isSnp(variant) || isDeletion(variant)) {
              val (_, allele, obs) = observed.head
              if (observed.count(_._2.nonEmpty) == 1 &&
                allele == variant.getAlternateAllele) {
                Some((variant, obs.duplicate(Some(false))))
              } else if (!obs.isRef ||
                observed.size != variant.getReferenceAllele.length) {
                Some((variant, obs.nullOut))
              } else {
                Some((variant, obs.invert))
              }
            } else if (isInsertion(variant)) {
              val insAllele = variant.getAlternateAllele.tail
              val insObserved = observed.filter(_._2 == insAllele)
              if (observed.size == 2 &&
                insObserved.size == 1) {
                Some((variant, insObserved.head._3.duplicate(Some(false))))
              } else if (observed.forall(_._3.isRef)) {
                Some((variant, observed.head._3.invert))
              } else {
                Some((variant, observed.head._3.nullOut))
              }
            } else {
              None
            }
          })

          obsMap
        }
      } catch {
        case t: Throwable => ProcessException.time {
          log.error("Processing read %s failed with exception %s. Skipping...".format(
            read.getReadName, t))
          Iterable.empty
        }
      }
    }
  }

  private def isSnp(v: Variant): Boolean = {
    v.getReferenceAllele.length == 1 &&
      v.getAlternateAllele.length == 1
  }

  private def isDeletion(v: Variant): Boolean = {
    v.getReferenceAllele.length > 1 &&
      v.getAlternateAllele.length == 1
  }

  private def isInsertion(v: Variant): Boolean = {
    v.getReferenceAllele.length == 1 &&
      v.getAlternateAllele.length > 1
  }

  /**
   * Scores the putative variants covered by a read against a single read.
   *
   * @param rdd An RDD containing the product of joining variants against reads
   *   and then grouping by the reads.
   * @param ploidy The assumed copy number at each site.
   * @return Returns an RDD of (Variant, Observation) pairs.
   */
  private[genotyping] def readsToObservations(
    rdd: RDD[(AlignmentRecord, Iterable[Variant])],
    ploidy: Int): RDD[(Variant, Observation)] = ObserveReads.time {

    rdd.flatMap(readToObservations(_, ploidy))
      .reduceByKey(_.merge(_))
  }

  /**
   * Turns observations of variants into genotype calls.
   *
   * @param rdd RDD of (variant, observation) pairs to transform.
   * @return Returns an RDD of genotype calls.
   */
  private[genotyping] def observationsToGenotypes(
    rdd: RDD[(Variant, Observation)],
    sample: String): RDD[Genotype] = EmitGenotypes.time {

    rdd.map(observationToGenotype(_, sample))
  }

  private val TEN_DIV_LOG10 = 10.0 / mathLog(10.0)

  /**
   * @param logLikelihood The log likelihoods of an observed call.
   * @return Returns a tuple containing the (genotype state, quality).
   */
  private[genotyping] def genotypeStateAndQuality(
    logLikelihoods: Array[Double]): (Int, Double) = {

    // even for copy number 1, we have at least 2 likelihoods
    assert(logLikelihoods.length >= 2)

    @tailrec def getMax(iter: Iterator[Double],
                        idx: Int,
                        maxIdx: Int,
                        max: Double,
                        second: Double): (Int, Double, Double) = {
      if (!iter.hasNext) {
        (maxIdx, max, second)
      } else {

        // do we have a new max? if so, replace max and shift to second
        val currValue = iter.next
        val (nextMaxIdx, nextMax, nextSecond) = if (currValue > max) {
          (idx, currValue, max)
        } else if (currValue > second) {
          (maxIdx, max, currValue)
        } else {
          (maxIdx, max, second)
        }

        getMax(iter, idx + 1, nextMaxIdx, nextMax, nextSecond)
      }
    }

    // which of first two genotype likelihoods is higher?
    val (startMaxIdx, startMax, startSecond) =
      if (logLikelihoods(0) >= logLikelihoods(1)) {
        (0, logLikelihoods(0), logLikelihoods(1))
      } else {
        (1, logLikelihoods(1), logLikelihoods(0))
      }

    // get the max state and top two likelihoods
    val (state, maxLikelihood, secondLikelihood) =
      getMax(logLikelihoods.toIterator.drop(2),
        2,
        startMaxIdx,
        startMax,
        startSecond)

    // phred quality is 10 * (max - second) / log10
    val quality = TEN_DIV_LOG10 * (maxLikelihood - secondLikelihood)

    (state, quality)
  }

  /**
   * Turns a single variant observation into a genotype call.
   *
   * @param variant Tuple of (variant, observation) to transform into variant
   *   calls at the site.
   * @return Returns a called genotype.
   */
  private[genotyping] def observationToGenotype(
    variant: (Variant, Observation),
    sample: String): Genotype = {

    // unpack tuple
    val (v, obs) = variant

    // merge the log likelihoods
    val coverage = obs.alleleCoverage + obs.otherCoverage
    val logLikelihoods = obs.alleleLogLikelihoods

    // get the genotype state and quality
    val (genotypeState, qual) = genotypeStateAndQuality(logLikelihoods)

    // build the genotype call array
    val alleles = Seq.fill(genotypeState)({
      GenotypeAllele.ALT
    }) ++ Seq.fill(obs.copyNumber - genotypeState)({
      GenotypeAllele.REF
    })

    // set up strand bias seq and calculate p value
    val sbComponents = Seq(obs.otherForwardStrand,
      obs.otherCoverage - obs.otherForwardStrand,
      obs.alleleForwardStrand,
      obs.alleleCoverage - obs.alleleForwardStrand)
    val sbPValue = fisher(sbComponents(0), sbComponents(1),
      sbComponents(2), sbComponents(3))

    // add variant calling annotations
    val vcAnnotations = VariantCallingAnnotations.newBuilder
      .setRmsMapQ(sqrt(obs.squareMapQ / coverage.toDouble).toFloat)
      .setFisherStrandBiasPValue(sbPValue)
      .build

    Genotype.newBuilder()
      .setVariant(v)
      .setVariantCallingAnnotations(vcAnnotations)
      .setStart(v.getStart)
      .setEnd(v.getEnd)
      .setContigName(v.getContigName)
      .setSampleId(sample)
      .setStrandBiasComponents(sbComponents
        .map(i => i: java.lang.Integer))
      .setReadDepth(obs.totalCoverage)
      .setReferenceReadDepth(obs.otherCoverage)
      .setAlternateReadDepth(obs.alleleCoverage)
      .setGenotypeLikelihoods(logLikelihoods.map(d => d.toFloat: java.lang.Float)
        .toSeq)
      .setAlleles(alleles)
      .setGenotypeQuality(qual.toInt)
      .build
  }

  /**
   * @param n The number to compute the factorial of.
   * @param factorial The running factorial value. Do not provide!
   * @return the factorial in log space.
   */
  private[genotyping] def logFactorial(n: Int, factorial: Double = 0.0): Double = {
    assert(n >= 0)
    if (n <= 1) {
      factorial
    } else {
      logFactorial(n - 1, mathLog(n) + factorial)
    }
  }

  /**
   * Performs Fisher's exact test on strand bias observations.
   *
   * @param otherFwd The number of reads that don't support this allele that are
   *   mapped on the forward strand.
   * @param otherRef The number of reads that don't support this allele that are
   *   mapped on the reverse strand.
   * @param alleleFwd The number of reads that support this allele that are
   *   mapped on the forward strand.
   * @param alleleRef The number of reads that support this allele that are
   *   mapped on the reverse strand.
   * @return Returns the Phred scaled P-value for whether there is a significant
   *   difference between the strand biases for the two alleles.
   */
  private[genotyping] def fisher(otherFwd: Int,
                                 otherRev: Int,
                                 alleleFwd: Int,
                                 alleleRev: Int): Float = {

    val numerator = (logFactorial(otherFwd + otherRev) +
      logFactorial(otherFwd + alleleFwd) +
      logFactorial(alleleFwd + alleleRev) +
      logFactorial(otherRev + alleleRev))
    val denominator = (logFactorial(otherFwd) +
      logFactorial(otherRev) +
      logFactorial(alleleFwd) +
      logFactorial(alleleRev) +
      logFactorial(otherFwd +
        otherRev +
        alleleFwd +
        alleleRev))

    LogPhred.logErrorToPhred(numerator - denominator).toFloat
  }
}
