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

import org.bdgenomics.formats.avro.{
  ADAMContig,
  ADAMGenotype,
  ADAMGenotypeAllele,
  ADAMRecord,
  ADAMVariant
}
import scala.collection.immutable.{ SortedSet, TreeSet }
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
import scala.collection.mutable.HashSet
import scala.math._

object ReadCallAssemblyPhaser extends VariantCallCompanion {

  val callName = "AssemblyPhaser"
  val debug = false

  def apply(stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration,
            partitions: PartitionSet): ReadCallAssemblyPhaser = {

    // get config values
    val flankLength = config.getInt("flankLength", 40)
    val kmerLength = config.getInt("kmerLength", 20)
    val trimSpurs = config.getBoolean("trimSpurs", true)
    val trimThreshold: Option[Double] = if (config.containsKey("lowCoverageTrimmingThreshold")) {
      Some(config.getDouble("lowCoverageTrimmingThreshold"))
    } else {
      None
    }
    val maxEntries = config.getInt("maxEntries", 5)

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

    new ReadCallAssemblyPhaser(partitions,
      kmerLength,
      flankLength,
      maxEntries,
      trimSpurs,
      trimThreshold,
      haplotypeAlignerConfig,
      readAlignerConfig)
  }

}

/**
 * Phase (diploid) haplotypes with kmer assembly on active regions.
 */
class ReadCallAssemblyPhaser(val _partitions: PartitionSet,
                             val kmerLen: Int = 20,
                             val _flankLength: Int = 40,
                             val maxEntries: Int = 5,
                             val trimSpurs: Boolean = true,
                             val lowCoverageTrimmingThreshold: Option[Double] = None,
                             val _haplotypeAlignerConfig: TransitionMatrixConfiguration = TransitionMatrixConfiguration(),
                             val _readAlignerConfig: TransitionMatrixConfiguration = TransitionMatrixConfiguration(),
                             val weightThreshold: Double = 0.1,
                             val minAlts: Int = 2) extends ReadCallHaplotypes(
  _partitions,
  _flankLength,
  _haplotypeAlignerConfig,
  _readAlignerConfig) {

  assert(minAlts >= 1, "Cannot set min alt count to less than 1.")

  override val companion: VariantCallCompanion = ReadCallAssemblyPhaser

  /**
   * Performs assembly over a region. See:
   *
   * C.A. Albers, G. Lunter, D.G. MacArthur, G. McVean, W.H. Ouwehand, R. Durbin.
   * "Dindel: Accurate indel calls from short-read data." Genome Research 21 (2011).
   *
   * @param region Sequence of reads spanning the region.
   * @param reference String representing reference over the region.
   * @return Haplotype strings corresponding to the region.
   */
  override def generateHaplotypes(region: Seq[RichADAMRecord], reference: String): Seq[String] = {
    val readLen = region(0).getSequence.length
    val regionLen = reference.length

    val tg = System.nanoTime
    var kmerGraph = KmerGraph(kmerLen,
      readLen,
      regionLen,
      reference,
      region,
      flankLength,
      maxEntries,
      trimSpurs,
      lowCoverageTrimmingThreshold)

    val paths: scala.collection.mutable.SortedSet[KmerPath] = try {
      kmerGraph.allPaths
    } catch {
      case soe: java.lang.StackOverflowError => {
        log.error("Caught stack overflow error when assembling region on " +
          region.head.getContig.getContigName + " from " +
          region.head.getStart + " to " + region.last.end.get +
          ". Will only compute reference likelihoods here.")

        scala.collection.mutable.SortedSet()
      }
    }

    @tailrec def dropAlts(neededWeight: Double,
                          minHaplotypes: Int,
                          haplotypes: Seq[KmerPath],
                          pickedHaplotypes: Seq[KmerPath] = Seq()): Seq[KmerPath] = {
      // return if we're down to the minimum number of haplotypes, or if we would take
      // the current head haplotype
      if (haplotypes.length == 0 ||
        (pickedHaplotypes.length >= minHaplotypes &&
          haplotypes.head.meanWeight < neededWeight)) {
        pickedHaplotypes
      } else {
        dropAlts(neededWeight,
          minHaplotypes,
          haplotypes.drop(1),
          pickedHaplotypes ++ haplotypes.take(1))
      }
    }

    val maxWeight = paths.map(_.meanWeight).fold(0.0)(_ max _)
    val keepWeight = (1.0 - weightThreshold) * maxWeight

    dropAlts(keepWeight, minAlts, paths.toSeq.sortBy(_.meanWeight).reverse)
      .map(_.haplotypeString)
  }
}
