/*
 * Copyright (c) 2014. Regents of the University of California
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

import net.sf.samtools.{ CigarElement, CigarOperator }
import org.bdgenomics.adam.avro.{
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
import scala.collection.immutable.{ SortedSet, TreeSet }
import scala.collection.mutable.HashSet
import scala.math._

object ReadCallHaplotypesFromReads extends VariantCallCompanion {

  val callName = "CallHaplotypesFromReads"
  val debug = false

  def apply(stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration,
            partitions: PartitionSet): ReadCallHaplotypes = {

    // get config values
    val flankLength = config.getInt("flankLength", 40)

    // get aligner configs
    val haplotypeAlignerConfig = try {
      TransitionMatrixConfiguration(config.configurationAt("haplotypeAligner"))
    } catch {
      case _: Throwable => TransitionMatrixConfiguration()
    }
    val readAlignerConfig = try {
      TransitionMatrixConfiguration(config.configurationAt("readAligner"))
    } catch {
      case _: Throwable => TransitionMatrixConfiguration()
    }

    new ReadCallHaplotypesFromReads(partitions,
      flankLength,
      haplotypeAlignerConfig,
      readAlignerConfig)
  }
}

/**
 * Phase haplotypes with haplotypes extracted from reads.
 */
class ReadCallHaplotypesFromReads(
  val _partitions: PartitionSet,
  val _flankLength: Int = 40,
  val _haplotypeAlignerConfig: TransitionMatrixConfiguration = TransitionMatrixConfiguration(),
  val _readAlignerConfig: TransitionMatrixConfiguration = TransitionMatrixConfiguration()) extends ReadCallHaplotypes(
  _partitions,
  _flankLength,
  _haplotypeAlignerConfig,
  _readAlignerConfig) {

  val companion: VariantCallCompanion = ReadCallHaplotypesFromReads

  /**
   * Generate haplotypes to score. Here, we generate haplotypes by inserting read sequences
   * directly into the reference.
   *
   * @param region Reads to generate haplotypes from.
   * @param reference
   */
  def generateHaplotypes(region: Seq[RichADAMRecord], reference: String): Seq[String] = {
    // get start and end of region
    val start = region.map(_.getStart).min
    val end = region.flatMap(_.end).max

    // insert reads into the reference
    region.map(insertReadIntoReference(_, reference, start, end))
      .distinct
  }
}
