/*
 * Copyright (c) 2013-2014. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

import org.bdgenomics.adam.avro.{
  ADAMContig,
  ADAMGenotype,
  ADAMGenotypeAllele,
  ADAMRecord,
  ADAMVariant
}
import scala.collection.immutable.{ SortedSet, TreeSet }
import org.bdgenomics.adam.models.ADAMVariantContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichADAMRecord
import org.bdgenomics.adam.rich.RichADAMRecord._
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.algorithms.debrujin._
import org.bdgenomics.avocado.algorithms.hmm._
import org.bdgenomics.avocado.calls.VariantCallCompanion
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.{ RDD }
import scala.collection.mutable.HashSet
import scala.math._

object VariantType extends scala.Enumeration {
  type VariantType = Value
  val SNP, MNP, Insertion, Deletion = Value
}

object ReadCallAssemblyPhaser extends VariantCallCompanion {

  val callName = "AssemblyPhaser"
  val debug = false

  def apply(stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): ReadCallAssemblyPhaser = {

    new ReadCallAssemblyPhaser()
  }

}

/**
 * Phase (diploid) haplotypes with kmer assembly on active regions.
 */
class ReadCallAssemblyPhaser(val kmerLen: Int = 20,
                             val regionWindow: Int = 200,
                             val flankLength: Int = 40) extends ReadCall {

  val companion = ReadCallAssemblyPhaser

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
    // TODO(peter, 12/5) currently, get the reference subsequence from the
    // MD tags of the ADAM records. Not entirely correct, because a region may
    // not be completely covered by reads, in which case the MD tags are
    // insufficient, so we ultimately want to resolve using the ref itself,
    // and not just the reads.
    val posRefs =
      region.map(read => (read.getStart, read.mdTag.map(_.getReference(read))))
        .sortBy(_._1)
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
      if (readRef.isDefined) {
        val readReference = readRef.get

        val relPos = position - startPos
        val offset = reference.length - relPos.toInt
        if (offset >= 0 && offset < readReference.length) {
          try {
            reference += readReference.substring(offset)
          } catch {
            case (e: StringIndexOutOfBoundsException) => {
              log.warn("String index out of bounds at: " + reference + ", " + readReference + ", " + offset)
            }
          }
        } else if (offset < 0) {
          return ""
        }
      } else {
        log.warn("No reference discovered for read due to missing MDTag or unmapped read")
      }
    }
    reference
  }

  /**
   * Checks to see if region is active.
   *
   * @param region Sequence of reads over which to test.
   * @param ref Reference sequence over which to test.
   * @return True if region is active.
   */
  def isRegionActive(region: Seq[RichADAMRecord], ref: String): Boolean = {
    // TODO(peter, 12/6) a very naive active region criterion. Upgrade asap!
    val activeLikelihoodThresh = -2.0
    var refHaplotype = new Haplotype(ref, region)
    val readsLikelihood = refHaplotype.readsLikelihood
    readsLikelihood < activeLikelihoodThresh
  }

  /**
   * Performs assembly over a region.
   *
   * @param region Sequence of reads spanning the region.
   * @param reference String representing reference over the region.
   * @return Kmer graph corresponding to region.
   */

  def assemble(region: Seq[RichADAMRecord], reference: String, removeSpurs: Boolean = false): KmerGraph = {
    val readLen = region(0).getSequence.length
    val regionLen = min(regionWindow + readLen - 1, reference.length)
    var kmerGraph = KmerGraph(kmerLen, readLen, regionLen, reference, region, flankLength, removeSpurs)
    kmerGraph
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
      print("not calling")
      List[ADAMGenotype]()
    }
  }

  /**
   * Phasing the assembled haplotypes to call variants. See:
   *
   * C.A. Albers, G. Lunter, D.G. MacArthur, G. McVean, W.H. Ouwehand, R. Durbin.
   * "Dindel: Accurate indel calls from short-read data." Genome Research 21 (2011).
   *
   * @param region Sequence of reads covering region.
   * @param kmerGraph Graph of kmers to use for calling.
   * @param reference String for reference in the region.
   * @return List of variant contexts found in the region.
   */
  def phaseAssembly(region: Seq[RichADAMRecord],
                    kmerGraph: KmerGraph,
                    reference: String,
                    maxHaplotypes: Int = 16): List[ADAMVariantContext] = {

    val aligner = new HMMAligner
    val refHaplotype = new Haplotype(reference, region, aligner, reference)

    // Score all haplotypes against the reads.
    val orderedHaplotypes = SortedSet[Haplotype](kmerGraph.allPaths.map(path =>
      new Haplotype(path.haplotypeString, region, aligner, reference)).toSeq: _*)(HaplotypeOrdering.reverse)

    // Pick the top X-1 haplotypes and the reference haplotype.
    val bestHaplotypes = refHaplotype :: orderedHaplotypes.take(maxHaplotypes - 1).toList

    // Score the haplotypes pairwise inclusively.
    val orderedHaplotypePairBuilder = TreeSet.newBuilder[HaplotypePair](HaplotypePairOrdering.reverse)
    for (i <- 0 until bestHaplotypes.size) {
      for (j <- i until bestHaplotypes.size) {
        var pair = new HaplotypePair(bestHaplotypes(i), bestHaplotypes(j), aligner)
        orderedHaplotypePairBuilder += pair
      }
    }

    val orderedHaplotypePairs = orderedHaplotypePairBuilder.result

    if (ReadCallAssemblyPhaser.debug) {
      println("After scoring, have:")
      orderedHaplotypePairs.foreach(println)
    }

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
    if (calledHaplotypePair.isDefined) {
      var variants = List[ADAMGenotype]()
      val sortedRegion = region.sortBy(_.getStart)
      val firstRead = sortedRegion(0)
      val refPos = firstRead.getStart
      val refName = firstRead.getContig.getContigName.toString
      val sampleName = firstRead.getRecordGroupSample.toString
      var calledHaplotypes = new HashSet[Haplotype]
      val calledHaplotype1 = calledHaplotypePair.get.haplotype1
      val calledHaplotype2 = calledHaplotypePair.get.haplotype2
      if (ReadCallAssemblyPhaser.debug) {
        println("Called:")
        println(calledHaplotype1 + ", " + calledHaplotype1.hasVariants + ", " + calledHaplotype1.referenceAlignment.alignment)
        println(calledHaplotype2 + ", " + calledHaplotype2.hasVariants + ", " + calledHaplotype2.referenceAlignment.alignment)
      }
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
              val variant = emitVariantCall(variantType, variantLength,
                variantOffset, refOffset,
                haplotype.sequence, refHaplotype.sequence,
                heterozygousRef,
                heterozygousNonref,
                variantPhred,
                refPos,
                sampleName, refName)
              variants = variants ::: variant
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
      genotypesToVariantContext(variants)
    } else {
      List()
    }
  }

  /**
   * Method signature for variant calling operation.
   *
   * @param[in] pileupGroups An RDD containing reads.
   * @return An RDD containing called variants.
   */
  override def call(reads: RDD[ADAMRecord]): RDD[ADAMVariantContext] = {
    log.info("Grouping reads into active regions.")
    val richReads = reads.map(RichADAMRecord(_))
    val regions = richReads.groupBy(r => r.getStart / regionWindow)

    val regionsAndReference = regions.map(kv => (getReference(kv._2), kv._2))
    val activeRegions = regionsAndReference.filter(kv => isRegionActive(kv._2, kv._1))

    log.info("Calling variants with local assembly.")

    activeRegions.flatMap(x => {
      val ref = x._1
      val region = x._2
      if (ref.length > 0 && region.length > 0) {
        val kmerGraph = assemble(region, ref)
        phaseAssembly(region, kmerGraph, ref)
      } else {
        List()
      }
    })
  }

  override def isCallable(): Boolean = true
}
