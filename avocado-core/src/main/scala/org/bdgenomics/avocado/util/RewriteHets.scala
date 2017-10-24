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
package org.bdgenomics.avocado.util

import org.bdgenomics.adam.rdd.variant.GenotypeRDD
import org.bdgenomics.formats.avro.{ Genotype, GenotypeAllele }
import scala.collection.JavaConversions._

private[avocado] trait RewriteHetsArgs extends Serializable {

  /**
   * The maximum allele fraction for the alternate allele in a het SNP call.
   *
   * Set to a negative value to omit.
   */
  var maxHetSnpAltAllelicFraction: Float

  /**
   * The maximum allele fraction for the alternate allele in a het SNP call.
   *
   * Set to a negative value to omit.
   */
  var maxHetIndelAltAllelicFraction: Float

  /**
   * If true, does not attempt to rewrite het SNPs.
   */
  var disableHetSnpRewriting: Boolean

  /**
   * If true, does not attempt to rewrite het INDELs.
   */
  var disableHetIndelRewriting: Boolean
}

/**
 * Rewrites high allelic fraction het genotypes as homozygous alternate calls.
 */
object RewriteHets extends Serializable {

  /**
   * Identifies high allelic fraction het calls in an RDD of genotypes and
   * rewrites them as homozygous alt calls.
   *
   * @param rdd The RDD of genotypes to filter.
   * @param args The arguments to configure the rewriter.
   * @return Returns a new RDD of genotypes.
   */
  def apply(rdd: GenotypeRDD,
            args: RewriteHetsArgs): GenotypeRDD = {

    val maxSnpAllelicFraction = args.maxHetSnpAltAllelicFraction
    val maxIndelAllelicFraction = args.maxHetIndelAltAllelicFraction
    val rewriteHetSnps = !args.disableHetSnpRewriting
    val rewriteHetIndels = !args.disableHetIndelRewriting

    if (rewriteHetSnps || rewriteHetIndels) {
      rdd.transform(gtRdd => gtRdd.map(processGenotype(_,
        maxSnpAllelicFraction,
        maxIndelAllelicFraction,
        rewriteHetSnps,
        rewriteHetIndels)))
    } else {
      rdd
    }
  }

  /**
   * Examines a single genotype to see if it should be rewritten.
   *
   * @param gt The genotype to examine.
   * @param maxSnpAllelicFraction The threshold for considering a het SNP call
   *   to be a miscalled homozygous SNP.
   * @param maxIndelAllelicFraction The threshold for considering a het INDEL
   *   call to be a miscalled homozygous INDEL
   * @param rewriteHetSnps If false, disables SNP checking.
   * @param rewriteHetIndels If false, disables INDEL checking.
   * @return Returns true if the genotype should be rewritten..
   */
  private[util] def shouldRewrite(gt: Genotype,
                                  maxSnpAllelicFraction: Float,
                                  maxIndelAllelicFraction: Float,
                                  rewriteHetSnps: Boolean,
                                  rewriteHetIndels: Boolean): Boolean = {
    if (gt.getVariant.getAlternateAllele == null) {
      false
    } else {

      val isSnp = ((gt.getVariant.getReferenceAllele.length == 1) &&
        (gt.getVariant.getAlternateAllele.length == 1))
      val numAlts = gt.getAlleles.count(_ == GenotypeAllele.ALT)
      val isHet = (numAlts != 0) && (numAlts != gt.getAlleles.length)

      def checkAf(af: Float): Boolean = {
        (Option(gt.getReadDepth), Option(gt.getAlternateReadDepth)) match {
          case (Some(dp), Some(alt)) => (alt.toFloat / dp.toFloat) >= af
          case _                     => false
        }
      }

      if (rewriteHetSnps && isSnp && isHet) {
        checkAf(maxSnpAllelicFraction)
      } else if (rewriteHetIndels && !isSnp && isHet) {
        checkAf(maxIndelAllelicFraction)
      } else {
        false
      }
    }
  }

  /**
   * Rewrites a het genotype as a hom alt call.
   *
   * @param gt The genotype to rewrite.
   * @return Returns the rewritten genotype.
   */
  private[util] def rewriteGenotype(gt: Genotype): Genotype = {
    val numAlleles = gt.getAlleles.length
    val newAlleles = List.fill(numAlleles) { GenotypeAllele.ALT }

    Genotype.newBuilder(gt)
      .setGenotypeQuality(null)
      .setAlleles(newAlleles)
      .build
  }

  /**
   * Processes a single genotype, and rewrites it if it appears to be a
   * miscalled hom alt.
   *
   * @param gt The genotype to examine.
   * @param maxSnpAllelicFraction The threshold for considering a het SNP call
   *   to be a miscalled homozygous SNP.
   * @param maxIndelAllelicFraction The threshold for considering a het INDEL
   *   call to be a miscalled homozygous INDEL
   * @param rewriteHetSnps If false, disables SNP checking.
   * @param rewriteHetIndels If false, disables INDEL checking.
   * @return Returns the rewritten genotype if the genotype should be rewritten,
   *   else returns the original genotype.
   */
  private[util] def processGenotype(gt: Genotype,
                                    maxSnpAllelicFraction: Float,
                                    maxIndelAllelicFraction: Float,
                                    rewriteHetSnps: Boolean,
                                    rewriteHetIndels: Boolean): Genotype = {
    if (shouldRewrite(gt,
      maxSnpAllelicFraction,
      maxIndelAllelicFraction,
      rewriteHetSnps,
      rewriteHetIndels)) {
      rewriteGenotype(gt)
    } else {
      gt
    }
  }
}
