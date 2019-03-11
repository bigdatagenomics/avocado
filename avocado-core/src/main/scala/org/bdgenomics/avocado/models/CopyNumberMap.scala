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
package org.bdgenomics.avocado.models

import org.apache.spark.SparkContext._
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.feature.FeatureDataset
import scala.math.{ max, min }

private[avocado] object CopyNumberMap extends Serializable {

  /**
   * Creates a copy number variant map containing no CNVs.
   *
   * @param basePloidy The ploidy of this sample.
   * @return Returns a copy number map without CNVs.
   */
  def empty(basePloidy: Int): CopyNumberMap = {

    CopyNumberMap(basePloidy,
      Map.empty[String, Seq[(ReferenceRegion, Int)]])
  }

  /**
   * Creates a copy number variant map from CNVs stored as features.
   *
   * @param basePloidy The ploidy of this sample.
   * @param features Dataset of features.
   * @return Returns a map containing copy number variants.
   */
  def apply(basePloidy: Int,
            features: FeatureDataset): CopyNumberMap = {

    val cnvMap = features.rdd
      .flatMap(f => f.getFeatureType match {
        case "DUP" => Some((ReferenceRegion.unstranded(f),
          basePloidy + 1))
        case "DEL" => Some((ReferenceRegion.unstranded(f),
          basePloidy - 1))
        case _ => None
      }).groupBy(_._1.referenceName)
      .map(kv => {
        val (referenceName, cnvs) = kv

        (referenceName, cnvs.toSeq.sortBy(_._1))
      }).collect.toMap

    CopyNumberMap(basePloidy,
      cnvMap)
  }
}

/**
 * An object that stores copy number variation.
 *
 * @param basePloidy The ploidy of this sample.
 * @param variantsByReference A map mapping reference names to the regions containing
 *   copy number variants. These regions are sorted per reference, and are in
 *   tuples with the observed copy number over that region.
 */
private[avocado] case class CopyNumberMap private (
    val basePloidy: Int,
    private[models] val variantsByReference: Map[String, Seq[(ReferenceRegion, Int)]]) {

  /**
   * @return The lowest copy number seen over all regions.
   */
  def minPloidy: Int = {
    variantsByReference.values
      .flatMap(s => s.map(_._2))
      .fold(basePloidy)(_ min _)
  }

  /**
   * @return The highest copy number seen over all regions.
   */
  def maxPloidy: Int = {
    variantsByReference.values
      .flatMap(s => s.map(_._2))
      .fold(basePloidy)(_ max _)
  }

  /**
   * Returns copy number variants that overlap a query region.
   *
   * @param rr The region to query for overlap.
   * @return Returns all variants that overlap this region.
   */
  def overlappingVariants(
    rr: ReferenceRegion): Iterable[(ReferenceRegion, Int)] = {

    variantsByReference.get(rr.referenceName)
      .fold(Iterable.empty[(ReferenceRegion, Int)])(i => {
        i.dropWhile(!_._1.overlaps(rr))
          .takeWhile(_._1.overlaps(rr))
      })
  }
}
