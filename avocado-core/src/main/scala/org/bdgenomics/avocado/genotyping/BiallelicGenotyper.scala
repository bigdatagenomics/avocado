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
import org.apache.spark.sql.functions._
import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.types.IntegerType
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
  HardLimiter,
  LogPhred,
  LogUtils,
  TreeRegionJoin
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
    val joinedRdd = JoinReadsAndVariants.time {
      TreeRegionJoin.joinAndGroupByRight(
        variants.rdd.keyBy(v => ReferenceRegion(v)),
        reads.rdd.flatMap(r => {
          ReferenceRegion.opt(r).map(rr => (rr, r))
        })).map(_.swap)
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
    copyNumber: Int): Iterable[(DiscoveredVariant, SummarizedObservation)] = ObserveRead.time {

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
            } else if (isInsertion(variant)) {
              val insAllele = variant.getAlternateAllele.tail
              val insObserved = observed.filter(_._2 == insAllele)
              if (observed.size == 2 &&
                insObserved.size == 1) {
                Some((DiscoveredVariant(variant),
                  insObserved.head._3))
              } else if (observed.forall(_._3.isRef)) {
                Some((DiscoveredVariant(variant), observed.head._3.asRef))
              } else {
                Some((DiscoveredVariant(variant), observed.head._3.nullOut))
              }
            } else {
              val (_, allele, obs) = observed.head
              if (observed.count(_._2.nonEmpty) == 1 &&
                allele == variant.getAlternateAllele) {
                if (isDeletion(variant)) {
                  Some((DiscoveredVariant(variant), obs.asAlt))
                } else {
                  Some((DiscoveredVariant(variant), obs))
                }
              } else if (!obs.isRef ||
                observed.size != variant.getReferenceAllele.length) {
                Some((DiscoveredVariant(variant), obs.nullOut))
              } else {
                Some((DiscoveredVariant(variant), obs.asRef))
              }
            }
          })

          // check for overlapping other-alts
          obsMap.map(p => {
            val (dv, observed) = p
            if (observed.isNonRef) {
              val overlappingObs = obsMap.filter(kv => kv._1.overlaps(dv))
                .map(_._2)

              if (overlappingObs.exists(o => !o.isRef && !o.isOther && !o.isNonRef)) {
                (dv, observed.otherAlt)
              } else {
                p
              }
            } else {
              p
            }
          })
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

    val observations = rdd.flatMap(readToObservations(_, ploidy))

    // convert to dataframe
    val sqlContext = SQLContext.getOrCreate(rdd.context)
    import sqlContext.implicits._
    val observationsDf = sqlContext.createDataFrame(observations)

    // flatten schema
    val flatFields = Seq(
      observationsDf("_1.contigName").as("contigName"),
      observationsDf("_1.start").as("start"),
      observationsDf("_1.referenceAllele").as("referenceAllele"),
      observationsDf("_1.alternateAllele").as("alternateAllele"),
      observationsDf("_2.isRef").as("isRef"),
      observationsDf("_2.forwardStrand").as("forwardStrand"),
      observationsDf("_2.optQuality").as("optQuality"),
      observationsDf("_2.mapQ").as("mapQ"),
      observationsDf("_2.isOther").as("isOther"),
      observationsDf("_2.isNonRef").as("isNonRef"))
    val flatObservationsDf = observationsDf.select(flatFields: _*)

    // create scored table and prepare for join
    val scoredDf = broadcast(ScoredObservation.createFlattenedScores(
      rdd.context, 93, 93, ploidy))

    // run the join
    val joinedObservationsDf = scoredDf.join(flatObservationsDf,
      Seq("isRef",
        "isOther",
        "isNonRef",
        "forwardStrand",
        "optQuality",
        "mapQ"))

    // run aggregation
    val aggCols = Seq(
      sum("alleleForwardStrand").as("alleleForwardStrand"),
      sum("otherForwardStrand").as("otherForwardStrand"),
      sum("squareMapQ").as("squareMapQ")) ++ (0 to ploidy).map(i => {
        val field = "referenceLogLikelihoods%d".format(i)
        sum(field).as(field)
      }).toSeq ++ (0 to ploidy).map(i => {
        val field = "alleleLogLikelihoods%d".format(i)
        sum(field).as(field)
      }).toSeq ++ (0 to ploidy).map(i => {
        val field = "otherLogLikelihoods%d".format(i)
        sum(field).as(field)
      }).toSeq ++ (0 to ploidy).map(i => {
        val field = "nonRefLogLikelihoods%d".format(i)
        sum(field).as(field)
      }).toSeq ++ Seq(
        sum("alleleCoverage").as("alleleCoverage"),
        sum("otherCoverage").as("otherCoverage"),
        sum("totalCoverage").as("totalCoverage"),
        first("isRef").as("isRef"))
    val aggregatedObservationsDf = joinedObservationsDf.groupBy("contigName",
      "start",
      "referenceAllele",
      "alternateAllele")
      .agg(aggCols.head, aggCols.tail: _*)

    // re-nest the output
    val firstField = struct(aggregatedObservationsDf("contigName"),
      aggregatedObservationsDf("start"),
      aggregatedObservationsDf("referenceAllele"),
      aggregatedObservationsDf("alternateAllele"))
    val secondFields = Seq(
      aggregatedObservationsDf("alleleForwardStrand")
        .cast(IntegerType)
        .as("alleleForwardStrand"),
      aggregatedObservationsDf("otherForwardStrand")
        .cast(IntegerType)
        .as("otherForwardStrand"),
      aggregatedObservationsDf("squareMapQ"), {
        val fields = (0 to ploidy).map(i => {
          aggregatedObservationsDf("referenceLogLikelihoods%d".format(i))
        })
        array(fields: _*).as("referenceLogLikelihoods")
      }, {
        val fields = (0 to ploidy).map(i => {
          aggregatedObservationsDf("alleleLogLikelihoods%d".format(i))
        })
        array(fields: _*).as("alleleLogLikelihoods")
      }, {
        val fields = (0 to ploidy).map(i => {
          aggregatedObservationsDf("otherLogLikelihoods%d".format(i))
        })
        array(fields: _*).as("otherLogLikelihoods")
      }, {
        val fields = (0 to ploidy).map(i => {
          aggregatedObservationsDf("nonRefLogLikelihoods%d".format(i))
        })
        array(fields: _*).as("nonRefLogLikelihoods")
      },
      aggregatedObservationsDf("alleleCoverage")
        .cast(IntegerType)
        .as("alleleCoverage"),
      aggregatedObservationsDf("otherCoverage")
        .cast(IntegerType)
        .as("otherCoverage"),
      aggregatedObservationsDf("totalCoverage")
        .cast(IntegerType)
        .as("totalCoverage"),
      aggregatedObservationsDf("isRef"))
    val aggregatedNestedDf = aggregatedObservationsDf.select(firstField,
      struct(secondFields: _*))

    // convert to a dataset, and finally an rdd
    aggregatedNestedDf.as[(DiscoveredVariant, Observation)]
      .rdd
      .map(p => {
        val (v, obs) = p
        (v.toVariant, obs)
      })
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
   * @return Returns a tuple containing the (number of alt, ref, and otheralt
   *   alleles, quality).
   */
  private[genotyping] def genotypeStateAndQuality(
    obs: Observation): (Int, Int, Int, Double) = {

    val states = obs.referenceLogLikelihoods.size
    val scoreArray = new Array[Double](states)

    def blend(array1: Array[Double],
              array2: Array[Double]) {
      var idx1 = 0
      var idx2 = states - 1

      while (idx1 < states) {
        scoreArray(idx1) = array1(idx1) + array2(idx2)
        idx1 += 1
        idx2 -= 1
      }
    }

    // score alt/ref
    blend(obs.alleleLogLikelihoods, obs.referenceLogLikelihoods)
    val (altRefState, altRefQual) = genotypeStateAndQuality(scoreArray)

    // score alt/otheralt
    blend(obs.alleleLogLikelihoods, obs.otherLogLikelihoods)
    val (altOtherState, altOtherQual) = genotypeStateAndQuality(scoreArray)

    // score otheralt/ref
    blend(obs.otherLogLikelihoods, obs.referenceLogLikelihoods)
    val (otherAltState, otherAltQual) = genotypeStateAndQuality(scoreArray)

    if (altRefQual >= altOtherQual && altRefQual >= otherAltQual) {
      (altRefState, states - altRefState - 1, 0, altRefQual)
    } else if (altOtherQual >= otherAltQual) {
      (altOtherState, 0, states - altOtherState - 1, altOtherQual)
    } else {
      (0, states - otherAltState - 1, otherAltState, otherAltQual)
    }
  }

  /**
   * @param logLikelihood The log likelihoods of an observed call.
   * @return Returns a tuple containing the (genotype state, quality).
   */
  private[genotyping] def genotypeStateAndQuality(
    logLikelihoods: Array[Double]): (Int, Double) = {

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

    // get the genotype state and quality
    val (alt, ref, other, qual) = genotypeStateAndQuality(obs)

    // build the genotype call array
    val alleles = Seq.fill(alt)({
      GenotypeAllele.ALT
    }) ++ Seq.fill(ref)({
      GenotypeAllele.REF
    }) ++ Seq.fill(other)({
      GenotypeAllele.OTHER_ALT
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

    // merge likelihoods
    val likelihoods = obs.alleleLogLikelihoods.length
    val gl = new Array[java.lang.Float](likelihoods)
    val ol = new Array[java.lang.Float](likelihoods)
    def mergeArrays(a1: Array[Double],
                    a2: Array[Double],
                    oa: Array[java.lang.Float]) {
      var idx1 = 0
      var idx2 = likelihoods - 1
      while (idx1 < likelihoods) {
        oa(idx1) = (a1(idx1) + a2(idx2)).toFloat
        idx1 += 1
        idx2 -= 1
      }
    }
    mergeArrays(obs.alleleLogLikelihoods, obs.referenceLogLikelihoods, gl)
    mergeArrays(obs.nonRefLogLikelihoods, obs.referenceLogLikelihoods, ol)

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
      .setGenotypeLikelihoods(gl.toSeq)
      .setNonReferenceLikelihoods(ol.toSeq)
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
