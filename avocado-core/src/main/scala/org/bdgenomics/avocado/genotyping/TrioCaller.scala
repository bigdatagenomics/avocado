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

import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variant.{ GenotypeRDD, VariantContextRDD }
import org.bdgenomics.formats.avro.{ Genotype, GenotypeAllele }
import scala.collection.JavaConversions._

/**
 * Variant caller for genotyping the offspring of two parents.
 */
object TrioCaller extends Serializable {

  private val PHASED_FROM_FIRST: java.util.List[GenotypeAllele] =
    Seq(GenotypeAllele.ALT, GenotypeAllele.REF)
  private val PHASED_FROM_SECOND: java.util.List[GenotypeAllele] =
    Seq(GenotypeAllele.REF, GenotypeAllele.ALT)

  /**
   * Extracts the sample ID from a set of reads.
   *
   * Requires that read groups are attached, and that all read groups are from
   * a single sample.
   *
   * @param rdd The reads to extract the sample ID from.
   * @return The sample ID.
   */
  def extractSampleId(rdd: AlignmentRecordRDD): String = {
    require(!rdd.recordGroups.isEmpty, "Record groups are empty.")
    val samples = rdd.recordGroups.recordGroups
      .map(rg => rg.sample)
      .distinct
    require(samples.size == 1,
      "Had multiple sample names (%s) attached to reads.".format(
        samples.mkString(", ")))

    samples.head
  }

  /**
   * Trio calls genotypes in a pedigree with two parents and one child.
   *
   * @param rdd RDD of base genotypes.
   * @param firstParentId The sample ID for the first parent.
   * @param secondParentId The sample ID for the second parent.
   * @param childId The sample ID for the child.
   * @return Returns the final genotypes.
   */
  def apply(rdd: GenotypeRDD,
            firstParentId: String,
            secondParentId: String,
            childId: String): GenotypeRDD = {

    apply(rdd.toVariantContexts,
      firstParentId,
      secondParentId,
      childId).toGenotypes
  }

  /**
   * Trio calls genotypes in a pedigree with two parents and one child.
   *
   * @param rdd RDD of base genotypes.
   * @param firstParentId The sample ID for the first parent.
   * @param secondParentId The sample ID for the second parent.
   * @param childId The sample ID for the child.
   * @return Returns the final genotypes.
   */
  private[genotyping] def apply(rdd: VariantContextRDD,
                                firstParentId: String,
                                secondParentId: String,
                                childId: String): VariantContextRDD = {
    rdd.transform(rdd => {
      rdd.filter(!filterRef(_))
        .map(processVariant(_, firstParentId, secondParentId, childId))
        .filter(!filterRef(_))
    })
  }

  /**
   * Filters out sites where all genotypes are ref calls.
   *
   * @param vc A variant site.
   * @return True if all calls are ref.
   */
  private[genotyping] def filterRef(vc: VariantContext): Boolean = {
    vc.genotypes.isEmpty || vc.genotypes.forall(gt => {
      gt.getAlleles.forall(allele => {
        allele != GenotypeAllele.ALT
      })
    })
  }

  /**
   * Updates genotype calls for a trio at a single site.
   *
   * @param vc The variant site.
   * @param firstParentId The sample ID for the first parent.
   * @param secondParentId The sample ID for the second parent.
   * @param childId The sample ID for the child.
   * @return Returns a new variant context with the final genotypes.
   */
  private[genotyping] def processVariant(vc: VariantContext,
                                         firstParentId: String,
                                         secondParentId: String,
                                         childId: String): VariantContext = {

    def makeNoCall(sampleId: String): Genotype = {
      Genotype.newBuilder
        .setContigName(vc.variant.variant.getContigName)
        .setStart(vc.variant.variant.getStart)
        .setEnd(vc.variant.variant.getEnd)
        .setVariant(vc.variant.variant)
        .setAlleles(Seq(GenotypeAllele.NO_CALL, GenotypeAllele.NO_CALL))
        .setSampleId(sampleId)
        .build
    }

    def replaceWithNoCall(gt: Genotype): Genotype = {
      Genotype.newBuilder(gt)
        .setAlleles(Seq(GenotypeAllele.NO_CALL, GenotypeAllele.NO_CALL))
        .setGenotypeQuality(null)
        .build
    }

    // do we have calls for each sample?
    val optFirstParentGt = vc.genotypes.find(_.getSampleId == firstParentId)
    val optSecondParentGt = vc.genotypes.find(_.getSampleId == secondParentId)
    val optChildGt = vc.genotypes.find(_.getSampleId == childId)

    (optFirstParentGt, optSecondParentGt, optChildGt) match {
      case (Some(firstParentGt), Some(secondParentGt), Some(childGt)) => {

        def checkAlleles(): Boolean = {
          def checkGt(gt: Genotype): Boolean = {
            (gt.getAlleles.size == 2 &&
              gt.getAlleles.count(_ == GenotypeAllele.NO_CALL) == 0 &&
              gt.getAlleles.count(_ == GenotypeAllele.OTHER_ALT) == 0)
          }
          checkGt(firstParentGt) && checkGt(secondParentGt) && checkGt(childGt)
        }

        if (checkAlleles()) {

          // compute genotype states (count of alts)
          val firstParentState = firstParentGt.getAlleles.count(_ == GenotypeAllele.ALT)
          val secondParentState = secondParentGt.getAlleles.count(_ == GenotypeAllele.ALT)
          val childState = childGt.getAlleles.count(_ == GenotypeAllele.ALT)

          val newChildGt = if (childState == 0) {
            if (firstParentState == 2 || secondParentState == 2) {
              replaceWithNoCall(childGt)
            } else {
              Genotype.newBuilder(childGt)
                .setPhased(true)
                .build
            }
          } else if (childState == 2) {
            if (firstParentState == 0 || secondParentState == 0) {
              replaceWithNoCall(childGt)
            } else {
              Genotype.newBuilder(childGt)
                .setPhased(true)
                .build
            }
          } else {
            if ((firstParentState == 0 && secondParentState == 0) ||
              (firstParentState == 2 && secondParentState == 2)) {
              replaceWithNoCall(childGt)
            } else if (firstParentState == 1 &&
              secondParentState == 1) {
              // cannot phase het
              childGt
            } else {
              // can phase het
              val newAlleles = if (firstParentState != 0) {
                PHASED_FROM_FIRST
              } else {
                PHASED_FROM_SECOND
              }
              Genotype.newBuilder(childGt)
                .setPhased(true)
                .setAlleles(newAlleles)
                .build
            }
          }

          VariantContext(vc.variant.variant,
            Iterable(firstParentGt, secondParentGt, newChildGt))
        } else {
          VariantContext(vc.variant.variant,
            Iterable(firstParentGt, secondParentGt, childGt))
        }
      }
      case _ => {
        // fill in missing samples and return
        VariantContext(vc.variant.variant,
          Iterable(
            optFirstParentGt.getOrElse(makeNoCall(firstParentId)),
            optSecondParentGt.getOrElse(makeNoCall(secondParentId)),
            optChildGt.getOrElse(makeNoCall(childId))))
      }
    }
  }
}
