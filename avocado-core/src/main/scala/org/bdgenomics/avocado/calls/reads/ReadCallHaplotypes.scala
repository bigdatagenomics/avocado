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
package org.bdgenomics.avocado.calls.reads

import net.sf.samtools.{ CigarElement, CigarOperator }
import org.bdgenomics.formats.avro.{
  ADAMContig,
  ADAMGenotype,
  ADAMGenotypeAllele,
  ADAMRecord,
  ADAMVariant
}
import org.bdgenomics.adam.models.{ ADAMVariantContext, ReferenceRegion }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichADAMRecord
import org.bdgenomics.adam.rich.RichADAMRecord._
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.algorithms.debrujin._
import org.bdgenomics.avocado.algorithms.hmm._
import org.bdgenomics.avocado.calls.VariantCallCompanion
import org.bdgenomics.avocado.partitioners.PartitionSet
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import scala.annotation.tailrec
import scala.collection.immutable.{ SortedSet, TreeSet }
import scala.collection.mutable.HashSet
import scala.math._

object VariantType extends scala.Enumeration {
  type VariantType = Value
  val SNP, MNP, Insertion, Deletion = Value
}

/**
 * Phase haplotypes with haplotypes extracted from reads.
 */
abstract class ReadCallHaplotypes(
    val partitions: PartitionSet,
    val flankLength: Int = 40,
    val haplotypeAlignerConfig: TransitionMatrixConfiguration = TransitionMatrixConfiguration(),
    val readAlignerConfig: TransitionMatrixConfiguration = TransitionMatrixConfiguration()) extends ReadCall {

  val companion: VariantCallCompanion

  /**
   * Gets the reference from a set of reads. Works _provided_ that region has non-zero coverage
   * across whole region.
   *
   * @param region Sequence containing reads that cover region.
   * @return String of reference bases covering the region.
   *
   * @see getReadReference
   */
  def getReference(region: Seq[RichADAMRecord]): String = {

    val ts = System.nanoTime

    // TODO(peter, 12/5) currently, get the reference subsequence from the
    // MD tags of the ADAM records. Not entirely correct, because a region may
    // not be completely covered by reads, in which case the MD tags are
    // insufficient, so we ultimately want to resolve using the ref itself,
    // and not just the reads.
    val posRefs = region.flatMap(read => {
      try {
        if (read.mdTag.isDefined) {
          Some((read.getStart, read.mdTag.get.getReference(read)))
        } else {
          log.warn("No reference recovered for read " + read.getReadName + " as MD tag was missing.")
          None
        }
      } catch {
        case _: Throwable => {
          log.warn("No reference recovered for read " + read.getReadName + " due to issue with MD tag.")
          None
        }
      }
    }).sortBy(_._1)
    val startPos = posRefs(0)._1
    var reference = ""
    for ((position, readRef) <- posRefs) {
      // Here's an explanatory picture:
      //
      // OK:   [-----ref-----)
      //             [---read---)
      //
      // Skip: [-----ref-----)
      //         [---read---)
      //
      // Bail: [-----ref-----)
      //                         [---read---)
      val readReference = readRef

      val relPos = position - startPos
      val offset = reference.length - relPos.toInt
      if (offset < 0) {
        reference += readReference
        log.warn("Have gap in reference at " + startPos)
      } else if (offset < readReference.length) {
        try {
          reference += readReference.substring(offset)
        } catch {
          case (e: StringIndexOutOfBoundsException) => {
            log.error("String index out of bounds at: " + reference + ", " + readReference + ", " + offset)
          }
        }
      }
    }

    reference
  }

  /**
   * Checks to see if region is active using an Illumina specific test. We mark
   * the region as active if the region shows indel evidence, or if the number
   * of mismatching bases in the region is higher than expected.
   *
   * @param region Sequence of reads over which to test.
   * @param ref Reference sequence over which to test.
   * @return True if region is active.
   */
  def isRegionActive(region: Seq[RichADAMRecord], ref: String): Boolean = {
    @tailrec def collectStats(hasIndel: Boolean,
                              mismatches: Int,
                              bases: Int,
                              reg: Iterator[RichADAMRecord]): (Boolean, Int, Int) = {
      if (hasIndel || !reg.hasNext) {
        (hasIndel, mismatches, bases)
      } else {
        val read = reg.next
        collectStats(read.record.getCigar.toString.contains('I') ||
          read.record.getCigar.toString.contains('D'),
          read.record.getMismatchingPositions.toString.count(_.isLetter) + mismatches,
          read.record.getSequence.length + bases,
          reg)
      }
    }

    val (hasIndel, mismatches, bases) = collectStats(false, 0, 0, region.toIterator)

    hasIndel || (mismatches.toDouble > bases.toDouble / 45.0)
  }

  /**
   * Generate haplotypes to score. Here, we generate haplotypes by inserting read sequences
   * directly into the reference.
   *
   * @param region Reads to generate haplotypes from.
   * @param reference
   */
  def generateHaplotypes(region: Seq[RichADAMRecord], reference: String): Seq[String]

  /**
   * Inserts the aligned sequence of a read into the reference.
   *
   * @param read Read to insert.
   * @param reference Reference string.
   * @param referenceStart The start position of this reference string.
   * @param referenceEnd The end position of this reference string.
   * @return Returns a haplotype.
   */
  def insertReadIntoReference(read: RichADAMRecord,
                              reference: String,
                              referenceStart: Long,
                              referenceEnd: Long): String = {
    // get read start and end
    val readStart = read.getStart
    val readEnd = read.end.get

    // find length of trimming necessary
    val cigar = read.samtoolsCigar
      .getCigarElements
      .filter(ce => ce.getOperator == CigarOperator.H) // hard clipping leads to read having already been trimmed
    val startClipping = cigar.takeWhile(ce => ce.getOperator == CigarOperator.S)
      .map(_.getLength)
      .fold(0)(_ + _)
    val endClipping = cigar.reverse
      .takeWhile(ce => ce.getOperator == CigarOperator.S)
      .map(_.getLength)
      .fold(0)(_ + _)

    // trim read for insertion into reference
    val trimmedRead = read.getSequence.toString.drop(startClipping).dropRight(endClipping)

    // get reference start and end
    val referencePrefix = reference.take((readStart - referenceStart).toInt)
    val referenceSuffix = reference.takeRight((referenceEnd - readEnd).toInt)

    // insert sequence into reference and return
    referencePrefix + trimmedRead + referenceSuffix
  }

  /**
   * Emits a variant call, if a sequence is either heterozygous or homozygous non-ref.
   *
   * @param varType Type of variant.
   * @param varLength Length of variant.
   * @param varOffset Offset of variant in assembled haplotype.
   * @param refOffset Offset of variant versus reference haplotype.
   * @param varSequence Assembled sequence for variant.
   * @param refSequence Reference sequence.
   * @param heterozygousRef True if variant is heterozygous with a reference allele.
   * @param heterozygousNonref True if variant is heterozygous with two non-reference alleles.
   * @param phred Phred scaled variant quality.
   * @param sampleName Name of sample.
   * @param refName Name of reference sequence.
   * @return List of genotypes.
   */
  def emitVariantCall(varType: VariantType.VariantType,
                      varLength: Int,
                      varOffset: Int,
                      refOffset: Int,
                      varSequence: String,
                      refSequence: String,
                      heterozygousRef: Boolean,
                      heterozygousNonref: Boolean,
                      phred: Int,
                      refPos: Long,
                      sampleName: String,
                      refName: String): List[ADAMGenotype] = {
    assert(!(heterozygousRef && heterozygousNonref))

    val refAllele = if (varType != VariantType.Insertion && varType != VariantType.Deletion) {
      refSequence.substring(refOffset, refOffset + varLength)
    } else if (varType == VariantType.Deletion) {
      refSequence.substring(refOffset - 1, refOffset + varLength)
    } else {
      refSequence.substring(refOffset, refOffset + 1)
    }
    val altAllele = if (varType != VariantType.Deletion && varType != VariantType.Insertion) {
      varSequence.substring(varOffset, varOffset + varLength)
    } else if (varType == VariantType.Insertion) {
      varSequence.substring(varOffset, varOffset + varLength + 1)
    } else {
      varSequence.substring(varOffset - 1, varOffset)
    }

    if (heterozygousRef) {
      val alleles = List(ADAMGenotypeAllele.Ref, ADAMGenotypeAllele.Alt)

      val contig = ADAMContig.newBuilder
        .setContigName(refName)
        .build
      val variant = ADAMVariant.newBuilder
        .setContig(contig)
        .setReferenceAllele(refAllele)
        .setVariantAllele(altAllele)
        .setPosition(refPos + refOffset)
        .build
      val genotype = ADAMGenotype.newBuilder()
        .setVariant(variant)
        .setSampleId(sampleName)
        .setGenotypeQuality(phred)
        .setExpectedAlleleDosage(1.0f)
        .setAlleles(alleles)
        .build()

      List(genotype)
    } else if (heterozygousNonref) {
      val alleles = List(ADAMGenotypeAllele.OtherAlt, ADAMGenotypeAllele.Alt)

      val contig = ADAMContig.newBuilder
        .setContigName(refName)
        .build
      val variant = ADAMVariant.newBuilder
        .setContig(contig)
        .setReferenceAllele(refAllele)
        .setVariantAllele(altAllele)
        .setPosition(refPos + refOffset)
        .build
      val genotype = ADAMGenotype.newBuilder()
        .setVariant(variant)
        .setSampleId(sampleName)
        .setGenotypeQuality(phred)
        .setExpectedAlleleDosage(1.0f)
        .setAlleles(alleles)
        .build()

      List(genotype)
    } else if (!heterozygousRef && !heterozygousNonref) {
      val alleles = List(ADAMGenotypeAllele.Alt, ADAMGenotypeAllele.Alt)

      val contig = ADAMContig.newBuilder
        .setContigName(refName)
        .build
      val variant = ADAMVariant.newBuilder
        .setContig(contig)
        .setReferenceAllele(refAllele)
        .setVariantAllele(altAllele)
        .setPosition(refPos + refOffset)
        .build
      val genotype = ADAMGenotype.newBuilder()
        .setVariant(variant)
        .setSampleId(sampleName)
        .setGenotypeQuality(phred)
        .setExpectedAlleleDosage(1.0f)
        .setAlleles(alleles)
        .build()

      List(genotype)
    } else {
      List[ADAMGenotype]()
    }
  }

  /**
   * Scoring the assembled haplotypes to call variants.
   *
   * @param region Sequence of reads covering region.
   * @param haplotypes Haplotype sequences to score.
   * @param reference String for reference in the region.
   * @param referenceRegion Describes the coordinates that this region overlaps.
   * @return List of variant contexts found in the region.
   */
  def scoreHaplotypes(region: Seq[RichADAMRecord],
                      haplotypes: Seq[String],
                      reference: String,
                      referenceRegion: Option[ReferenceRegion] = None,
                      maxHaplotypes: Int = 16): List[ADAMVariantContext] = {

    val ts = System.nanoTime

    // info for logging
    val (start, end) = referenceRegion.fold((region.map(_.getStart).min,
      region.flatMap(_.end).max))(rr => (rr.start, rr.end))
    val refName = region.head.getContig.getContigName

    val aligner = new HMMAligner(haplotypeAlignerConfig)
    val refHaplotype = new Haplotype(reference, region, aligner, reference, readAlignerConfig)

    // Score all haplotypes against the reads.
    val orderedHaplotypes = SortedSet[Haplotype](haplotypes.flatMap(path => {
      try {
        Some(new Haplotype(path, region, aligner, reference))
      } catch {
        case _: Throwable => None
      }
    }): _*)(HaplotypeOrdering.reverse)

    // print logging info
    if (orderedHaplotypes.size > 0) {
      val refLen = reference.length
      val hapNum = orderedHaplotypes.size
      val hapLens = orderedHaplotypes.map(_.sequence.length)
      val minLen = hapLens.min
      val maxLen = hapLens.max
      val avgLen = (hapLens.sum.toDouble / hapNum.toDouble).toInt
      val readCount = region.length
      log.info("In region with refName " + refName + " from " + start + " to " + end +
        ", have " + readCount + " reads and " + hapNum +
        " haplotypes with minimum length " + minLen +
        ", maximum length " + maxLen + ", average length " + avgLen +
        ", and reference length " + refLen + ".")
    } else {
      val readCount = region.length
      log.info("In region with refName " + refName + " from " + start + " to " + end +
        ", have " + readCount + " reads but no alt haplotypes.")

      // return early
      return List()
    }

    // Pick the top X-1 haplotypes and the reference haplotype.
    val bestHaplotypes = refHaplotype :: orderedHaplotypes.filter(_.hasVariants)
      .take(maxHaplotypes - 1)
      .toList

    // Score the haplotypes pairwise inclusively.
    val orderedHaplotypePairBuilder = TreeSet.newBuilder[HaplotypePair](HaplotypePairOrdering.reverse)
    for (i <- 0 until bestHaplotypes.size) {
      for (j <- i until bestHaplotypes.size) {
        var pair = new HaplotypePair(bestHaplotypes(i), bestHaplotypes(j), aligner)
        orderedHaplotypePairBuilder += pair
      }
    }

    val orderedHaplotypePairs = orderedHaplotypePairBuilder.result

    // Pick the best haplotype pairs with and without indels.
    val calledHaplotypePair = orderedHaplotypePairs.find(_.hasVariants)
    val uncalledHaplotypePair = orderedHaplotypePairs.find(!_.hasVariants).get

    // Compute the variant error probability and the equivalent phred score,
    // and use them for all called variants.
    val variantErrorProb = if (calledHaplotypePair.isDefined) {
      val calledProbability = pow(10.0, calledHaplotypePair.get.pairLikelihood)
      val uncalledProbability = pow(10.0, uncalledHaplotypePair.pairLikelihood)
      uncalledProbability / (calledProbability + uncalledProbability)
    } else {
      1.0
    }
    val variantPhred = PhredUtils.successProbabilityToPhred(variantErrorProb)

    // TODO(peter, 12/8) Call variants, _without_ incorporating phasing for now.
    val calls = if (calledHaplotypePair.isDefined) {
      var variants = List[ADAMVariantContext]()
      val sortedRegion = region.sortBy(_.getStart)
      val firstRead = sortedRegion(0)
      val refPos = firstRead.getStart
      val refName = firstRead.getContig.getContigName.toString
      val sampleName = firstRead.getRecordGroupSample.toString
      var calledHaplotypes = new HashSet[Haplotype]
      val calledHaplotype1 = calledHaplotypePair.get.haplotype1
      val calledHaplotype2 = calledHaplotypePair.get.haplotype2

      calledHaplotypes += calledHaplotypePair.get.haplotype1
      calledHaplotypes += calledHaplotypePair.get.haplotype2

      // FIXME this isn't quite right...
      val heterozygousRef = !calledHaplotypes.forall(_.hasVariants) && !calledHaplotypes.forall(!_.hasVariants)
      val heterozygousNonref = calledHaplotypes.filter(_.hasVariants).size > 1

      for (haplotype <- calledHaplotypes) {
        if (haplotype.hasVariants) {
          var variantOffset = 0
          var refOffset = 0
          for (tok <- haplotype.referenceAlignment.alignment) {
            val variantLength = tok._1
            val move = tok._2
            if (variantLength > 0 && (move == 'X' || move == 'I' || move == 'D')) {
              val variantType = move match {
                case 'X' => {
                  if (variantLength > 1) {
                    VariantType.MNP
                  } else {
                    VariantType.SNP
                  }
                }
                case 'I' => VariantType.Insertion
                case 'D' => VariantType.Deletion
              }
              try {
                val variant = emitVariantCall(variantType, variantLength,
                  variantOffset, refOffset,
                  haplotype.sequence, refHaplotype.sequence,
                  heterozygousRef,
                  heterozygousNonref,
                  variantPhred,
                  refPos,
                  sampleName, refName)
                variants = variants ::: genotypesToVariantContext(variant)
              } catch {
                case _: Throwable => {
                  log.warn("Variant call oddity experienced at " + (variantOffset + refPos) + " on " + refName)
                }
              }
            }
            if (move != 'D') {
              variantOffset += variantLength
            }
            if (move != 'I') {
              refOffset += variantLength
            }
          }
        }
      }
      val baseVariants = variants.length
      val filteredVariants = variants.filter(v => {
        val pos = v.position.pos
        pos >= start && pos < end
      })

      log.info("Called " + filteredVariants.length + " variants in region with refName " +
        refName + " from " + start + " to " + end + " (" + baseVariants + "with flank).")
      filteredVariants
    } else {
      log.info("Called no variants in region with refName " + refName + " from " +
        start + " to " + end + ".")
      List()
    }

    calls
  }

  /**
   * Method signature for variant calling operation.
   *
   * @param[in] pileupGroups An RDD containing reads.
   * @return An RDD containing called variants.
   */
  override def call(reads: RDD[ADAMRecord]): RDD[ADAMVariantContext] = {
    // broadcast partition set
    val bcastPset = reads.context.broadcast(partitions)

    // group reads
    log.info("Grouping reads into regions.")
    val richReads = reads.map(RichADAMRecord(_))
    val regions = richReads.flatMap(r => {
      val rr = ReferenceRegion(r.record).get
      val partitions = bcastPset.value.getPartition(rr, flankLength)
      partitions.map(p => (p, r))
    }).groupByKey()

    log.info("avocado: Have " + regions.count + " regions.")

    // generate region maps
    val regionsAndReference: RDD[(String, Seq[RichADAMRecord], Option[ReferenceRegion])] = regions.map(kv => {
      (getReference(kv._2.toSeq), kv._2.toSeq, bcastPset.value.getRegion(kv._1))
    })

    log.info("Calling variants with local assembly.")

    // TODO: add back active region filtering criteria
    regionsAndReference.flatMap(x => {
      val ref = x._1
      val region = x._2
      val refRegion = x._3
      if (ref.length > 0 && region.length > 0) {
        val haplotypes = generateHaplotypes(region, ref)
        scoreHaplotypes(region, haplotypes, ref, refRegion)
      } else if (ref.length == 0 && region.length > 0) {
        val start = region.map(_.getStart).min
        val end = region.flatMap(_.end).max
        val name = region.head.getContig.getContigName
        log.warn("Had null reference: " + name + " from " + start + " to " + end)
        List()
      } else {
        log.warn("Had null reads/reference. " + ref.length)
        List()
      }
    })
  }

  override def isCallable(): Boolean = true
}
