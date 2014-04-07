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
import org.bdgenomics.adam.models.ADAMVariantContext
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichADAMRecord
import org.bdgenomics.adam.rich.RichADAMRecord._
import org.bdgenomics.adam.util.{ MdTag, PhredUtils }
import org.bdgenomics.avocado.algorithms.debrujin._
import org.bdgenomics.avocado.algorithms.hmm._
import org.bdgenomics.avocado.calls.VariantCallCompanion
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import net.sf.samtools.{ Cigar, CigarOperator, CigarElement, TextCigarCodec }
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.{ SparkContext, Logging }
import org.apache.spark.rdd.{ RDD }
import scala.collection.JavaConversions._
import scala.collection.mutable.{ ArrayBuffer, Buffer, HashMap, HashSet, PriorityQueue, StringBuilder }
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
                             val regionWindow: Int = 200) extends ReadCall {

  val companion = ReadCallAssemblyPhaser

  /**
   * From a read, returns the reference sequence.
   *
   * @param read Read from which to return sequence.
   * @return String containing reference sequence over this read.
   *
   * @see https://github.com/bigdatagenomics/adam/blob/indel-realign/adam-commands/src/main/scala/edu/berkeley/cs/amplab/adam/util/MdTag.scala
   * @see getReference
   */
  def getReadReference(read: ADAMRecord): String = {
    val mdtag = MdTag(read.getMismatchingPositions.toString, read.getStart)

    val readSeq = RichADAMRecord(read).getSequence.toString
    val cigar = RichADAMRecord(read).samtoolsCigar

    var readPos = 0
    var refPos = 0
    var reference = ""

    val cigarEls: Buffer[CigarElement] = cigar.getCigarElements
    for (el <- cigarEls) {
      el.getOperator match {
        case CigarOperator.M => {
          for (i <- (0 until el.getLength)) {
            mdtag.mismatchedBase(refPos) match {
              case Some(b) => reference += b
              case None    => reference += readSeq(readPos)
            }
            readPos += 1
            refPos += 1
          }
        }
        case CigarOperator.D => {
          for (i <- (0 until el.getLength)) {
            mdtag.deletedBase(refPos) match {
              case Some(b) => reference += b
              case None    => {}
            }
            refPos += 1
          }
        }
        case CigarOperator.I => {
          readPos += el.getLength
        }
        case _ => {}
      }
    }

    reference
  }

  /**
   * Gets the reference from a set of reads. Works _provided_ that region has non-zero coverage
   * across whole region.
   *
   * @param region Sequence containing reads that cover region.
   * @return String of reference bases covering the region.
   *
   * @see getReadReference
   */
  def getReference(region: Seq[ADAMRecord]): String = {
    // TODO(peter, 12/5) currently, get the reference subsequence from the
    // MD tags of the ADAM records. Not entirely correct, because a region may
    // not be completely covered by reads, in which case the MD tags are
    // insufficient, so we ultimately want to resolve using the ref itself,
    // and not just the reads.
    val posRefs = (region.map(_.getStart), region.map(r => getReadReference(r)))
      .zipped.map((pos, ref) => (pos, ref))
      .sortBy(_._1)
    val startPos = posRefs(0)._1
    var reference = ""
    for ((pos, ref) <- posRefs) {
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
      val relPos = pos - startPos
      val offset = reference.length - relPos.toInt
      if (offset >= 0 && offset < ref.length) {
        try {
          reference += ref.substring(offset)
        }
        catch {
          case (e: StringIndexOutOfBoundsException) => {
            log.warn("String index out of bounds at: " + reference + ", " + ref + ", " + offset)
          }
        }
      }
      else if (offset < 0) {
        return ""
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
  def isRegionActive(region: Seq[ADAMRecord], ref: String): Boolean = {
    // TODO(peter, 12/6) a very naive active region criterion. Upgrade asap!
    val activeLikelihoodThresh = -2.0
    var refHaplotype = new Haplotype(ref)
    var hmm = new HMMAligner
    val readsLikelihood = refHaplotype.scoreReadsLikelihood(hmm, region)
    readsLikelihood < activeLikelihoodThresh
  }

  /**
   * Performs assembly over a region.
   *
   * @param region Sequence of reads spanning the region.
   * @param ref String representing reference over the region.
   * @return Kmer graph corresponding to region.
   */
  def assemble(region: Seq[ADAMRecord], ref: String): KmerGraph = {
    val readLen = region(0).getSequence.length
    val regionLen = min(regionWindow + readLen - 1, ref.length)
    var kmerGraph = new KmerGraph(kmerLen, readLen, regionLen, ref)
    kmerGraph.insertReads(region)
    kmerGraph.connectGraph
    //kmer_graph.removeSpurs // TODO(peter, 11/27) debug: not doing spur removal atm.
    kmerGraph.enumerateAllPaths
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
   * @param refId ID for reference.
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
                      refName: String,
                      refId: Int): List[ADAMGenotype] = {
    assert(!(heterozygousRef && heterozygousNonref))

    val refAllele = if (varType != VariantType.Insertion) {
      refSequence.substring(refOffset, refOffset + varLength)
    }
    else {
      ""
    }
    val altAllele = if (varType != VariantType.Deletion) {
      varSequence.substring(varOffset, varOffset + varLength)
    }
    else {
      ""
    }

    if (heterozygousRef) {
      val alleles = List(ADAMGenotypeAllele.Ref, ADAMGenotypeAllele.Alt)

      val contig = ADAMContig.newBuilder
        .setContigId(refId)
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
    }
    else if (!heterozygousRef && !heterozygousNonref) {
      val alleles = List(ADAMGenotypeAllele.Alt, ADAMGenotypeAllele.Alt)

      val contig = ADAMContig.newBuilder
        .setContigId(refId)
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
    }
    else {
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
   * @param ref String for reference in the region.
   * @return List of variant contexts found in the region.
   */
  def phaseAssembly(region: Seq[ADAMRecord],
                    kmerGraph: KmerGraph,
                    ref: String): List[ADAMVariantContext] = {
    var refHaplotype = new Haplotype(ref)

    // Score all haplotypes against the reads.
    var hmm = new HMMAligner
    refHaplotype.scoreReadsLikelihood(hmm, region)
    var orderedHaplotypes = new PriorityQueue[Haplotype]()(HaplotypeOrdering)
    for (path <- kmerGraph.getAllPaths) {
      val haplotype = new Haplotype(path.asHaplotypeString)
      haplotype.scoreReadsLikelihood(hmm, region)
      orderedHaplotypes.enqueue(haplotype)
    }

    // Pick the top X-1 haplotypes and the reference haplotype.
    val maxNumBestHaplotypes = 16
    val numBestHaplotypes = min(maxNumBestHaplotypes, orderedHaplotypes.length)
    var bestHaplotypes = new ArrayBuffer[Haplotype]
    bestHaplotypes += refHaplotype
    for (i <- 1 to numBestHaplotypes) {
      bestHaplotypes += orderedHaplotypes.dequeue
    }

    // Score the haplotypes pairwise inclusively.
    var refHaplotypePair: HaplotypePair = null
    var orderedHaplotypePairs = new PriorityQueue[HaplotypePair]()(HaplotypePairOrdering)
    for (i <- 0 until bestHaplotypes.length) {
      for (j <- i until bestHaplotypes.length) {
        var pair = new HaplotypePair(bestHaplotypes(i), bestHaplotypes(j))
        pair.scorePairLikelihood(hmm, region)
        orderedHaplotypePairs.enqueue(pair)
        if (i == 0 && j == 0) {
          refHaplotypePair = pair
        }
      }
    }

    if (ReadCallAssemblyPhaser.debug) {
      println("After scoring, have:")
      orderedHaplotypePairs.foreach(println)
    }

    // Pick the best haplotype pairs with and without indels.
    val (calledHaplotypePair, uncalledHaplotypePair) = {
      var calledRes: HaplotypePair = null
      var uncalledRes: HaplotypePair = null
      do {
        val res = orderedHaplotypePairs.dequeue
        res.alignToReference(hmm, refHaplotype) match {
          case true => {
            if (calledRes == null) {
              calledRes = res
            }
          }
          case false => {
            if (uncalledRes == null) {
              uncalledRes = res
            }
          }
        }
      } while ((calledRes == null || uncalledRes == null) && orderedHaplotypePairs.length > 0)
      // TODO(peter, 12/8) this ought to be a pathological bug if it ever
      // happens (i.e., the ref-ref pair counts as having variants).
      // On the other hand, there might not be any valid variant haplotypes.
      // (FIXME: Really I should be using Option[].)
      if (uncalledRes == null) {
        uncalledRes = refHaplotypePair
      }
      (calledRes, uncalledRes)
    }

    // Compute the variant error probability and the equivalent phred score,
    // and use them for all called variants.
    val variantErrorProb = if (calledHaplotypePair != null) {
      val calledProbability = pow(10.0, calledHaplotypePair.pairLikelihood)
      val uncalledProbability = pow(10.0, uncalledHaplotypePair.pairLikelihood)
      uncalledProbability / (calledProbability + uncalledProbability)
    }
    else {
      1.0
    }
    val variantPhred = PhredUtils.successProbabilityToPhred(variantErrorProb)

    // TODO(peter, 12/8) Call variants, _without_ incorporating phasing for now.
    if (calledHaplotypePair != null) {
      var variants = List[ADAMGenotype]()
      val sortedRegion = region.sortBy(_.getStart)
      val firstRead = sortedRegion(0)
      val refPos = firstRead.getStart
      val refName = firstRead.getReferenceName.toString
      val refId = firstRead.getReferenceId
      val sampleName = firstRead.getRecordGroupSample.toString
      var calledHaplotypes = new HashSet[Haplotype]
      val calledHaplotype1 = calledHaplotypePair.haplotype1
      val calledHaplotype2 = calledHaplotypePair.haplotype2
      if (ReadCallAssemblyPhaser.debug) {
        println("Called:")
        println(calledHaplotype1 + ", " + calledHaplotype1.hasVariants + ", " + calledHaplotype1.alignment)
        println(calledHaplotype2 + ", " + calledHaplotype2.hasVariants + ", " + calledHaplotype2.alignment)
      }
      calledHaplotypes += calledHaplotypePair.haplotype1
      calledHaplotypes += calledHaplotypePair.haplotype2
      // FIXME this isn't quite right...
      val heterozygousRef = !calledHaplotypes.forall(_.hasVariants) && !calledHaplotypes.forall(!_.hasVariants)
      val heterozygousNonref = calledHaplotypes.filter(_.hasVariants).size > 1
      for (haplotype <- calledHaplotypes) {
        if (haplotype.hasVariants) {
          var variantOffset = 0
          var refOffset = 0
          for (tok <- haplotype.alignment) {
            val variantLength = tok._1
            val move = tok._2
            if (variantLength > 0 && (move == 'X' || move == 'I' || move == 'D')) {
              val variantType = move match {
                case 'X' => {
                  if (variantLength > 1) {
                    VariantType.MNP
                  }
                  else {
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
                sampleName, refName, refId)
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
    }
    else {
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
    val activeRegions = reads.groupBy(r => r.getStart / regionWindow)
      .map(x => (getReference(x._2), x._2))
      .filter(x => isRegionActive(x._2, x._1))

    log.info("Calling variants with local assembly.")

    activeRegions.flatMap(x => {
      val ref = x._1
      val region = x._2
      if (ref.length > 0 && region.length > 0) {
        val kmerGraph = assemble(region, ref)
        phaseAssembly(region, kmerGraph, ref)
      }
      else {
        List()
      }
    })
  }

  override def isCallable(): Boolean = true
}
