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

import org.bdgenomics.formats.avro.{
  Genotype,
  Variant,
  VariantAnnotation
}
import scala.collection.JavaConversions._

/**
 * Helper object for creating variant summaries.
 */
private[genotyping] object VariantSummary {

  /**
   * Extracts a variant summary from a single genotype.
   *
   * @param g The genotype to extract variant level stats from.
   * @return Returns a variant summary from the variant calling annotations
   *   attached to this genotype.
   */
  def apply(g: Genotype): VariantSummary = {
    val readDepth = Option(g.getReadDepth)
      .map(i => i: Int)
    val referenceReadDepth = Option(g.getReferenceReadDepth)
      .map(i => i: Int)

    val (frd, rrd, frrd, rrrd) = if (g.getStrandBiasComponents.isEmpty) {
      (None, None, None, None)
    } else {
      val sbc: Seq[Int] = g.getStrandBiasComponents
        .map(i => i: Int)
      require(sbc.size == 4,
        "Genotype (%s) has an illegal number of strand bias components.".format(g))

      // strand bias component:
      // 0,1 --> ref
      // 2,3 --> alt
      // 0,2 --> fwd
      // 1,3 --> rev
      val rfs = sbc(0)
      val rrs = sbc(1)
      val afs = sbc(2)
      val ars = sbc(3)

      (Some(rfs + afs), Some(rrs + ars),
        Some(rfs), Some(rrs))
    }

    VariantSummary(readDepth, referenceReadDepth,
      frd, rrd,
      frrd, rrrd)
  }
}

/**
 * A roll-up of genotype stats to the variant level.
 *
 * @param readDepth The number of reads seen at this site.
 * @param referenceReadDepth The number of reads supporting the reference.
 * @param forwardReadDepth The number of reads mapped on the forward strand.
 * @param reverseReadDepth The number of reads mapped on the reverse strand.
 * @param forwardReferenceReadDepth The number of reads that support the
 *   reference that mapped on the forward strand.
 * @param reverseReferenceReadDepth The number of reads that support the
 *   reference that mapped on the reverse strand.
 */
private[genotyping] case class VariantSummary(
    readDepth: Option[Int],
    referenceReadDepth: Option[Int],
    forwardReadDepth: Option[Int],
    reverseReadDepth: Option[Int],
    forwardReferenceReadDepth: Option[Int],
    reverseReferenceReadDepth: Option[Int]) {

  private def mergeOptions(o1: Option[Int],
                           o2: Option[Int]): Option[Int] = {
    (o1, o2) match {
      case (Some(a), Some(b)) => Some(a + b)
      case (Some(a), None)    => Some(a)
      case (None, Some(b))    => Some(b)
      case (None, None)       => None
    }
  }

  /**
   * Merges two variant summaries together.
   *
   * @param vs The variant summary to merge into this summary.
   * @return Aggregates the statistics between the two summaries.
   */
  def merge(vs: VariantSummary): VariantSummary = {
    VariantSummary(
      mergeOptions(readDepth, vs.readDepth),
      mergeOptions(referenceReadDepth, vs.referenceReadDepth),
      mergeOptions(forwardReadDepth, vs.forwardReadDepth),
      mergeOptions(reverseReadDepth, vs.reverseReadDepth),
      mergeOptions(forwardReferenceReadDepth, vs.forwardReferenceReadDepth),
      mergeOptions(reverseReferenceReadDepth, vs.reverseReferenceReadDepth))
  }

  /**
   * Converts this summary into a variant annotation.
   *
   * @param variant The variant to base this annotation on. Used to get the
   *   current annotations on the variant.
   * @return Returns a new variant annotation object.
   */
  def toAnnotation(v: Variant): VariantAnnotation = {
    val vab = Option(v.getAnnotation).fold(VariantAnnotation.newBuilder)(va => {
      VariantAnnotation.newBuilder(va)
    })

    readDepth.map(i => i: java.lang.Integer)
      .foreach(vab.setReadDepth)
    referenceReadDepth.map(i => i: java.lang.Integer)
      .foreach(vab.setReferenceReadDepth)
    forwardReadDepth.map(i => i: java.lang.Integer)
      .foreach(vab.setForwardReadDepth)
    reverseReadDepth.map(i => i: java.lang.Integer)
      .foreach(vab.setReverseReadDepth)
    forwardReferenceReadDepth.map(i => i: java.lang.Integer)
      .foreach(vab.setReferenceForwardReadDepth)
    reverseReferenceReadDepth.map(i => i: java.lang.Integer)
      .foreach(vab.setReferenceReverseReadDepth)

    vab.build
  }
}
