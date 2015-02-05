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
package org.bdgenomics.avocado.discovery

import htsjdk.samtools.{ Cigar, CigarElement, CigarOperator, TextCigarCodec }
import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.SparkContext._
import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{
  ReferenceMapping,
  ReferencePosition,
  ReferenceRegion,
  SequenceDictionary
}
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.ShuffleRegionJoin
import org.bdgenomics.adam.rich.ReferenceMappingContext._
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.algorithms.debrujin.KmerGraph
import org.bdgenomics.avocado.models.{ Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.{ AlignmentRecord, NucleotideContigFragment }

object ReassemblyExplorer extends ExplorerCompanion with Serializable {

  val explorerName: String = "ReassemblyExplorer"
  @transient lazy private val CIGAR_CODEC: TextCigarCodec = TextCigarCodec.getSingleton

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Explorer = {
    new ReassemblyExplorer(config.getInt("kmerLength", 20),
      stats.reference,
      stats.sequenceDict,
      stats.contigLengths,
      config.getDouble("highCoverageThreshold", Double.PositiveInfinity),
      config.getDouble("lowCoverageThreshold", -1.0),
      config.getDouble("mismatchRateThreshold", 0.025),
      config.getDouble("clipRateThreshold", 0.05))
  }

  implicit object ContigReferenceMapping extends ReferenceMapping[(Long, NucleotideContigFragment)] with Serializable {
    override def getReferenceName(value: (Long, NucleotideContigFragment)): String = value._2.getContig.getContigName.toString
    override def getReferenceRegion(value: (Long, NucleotideContigFragment)): ReferenceRegion = ReferenceRegion(value._2).get
  }

  private[discovery] def calculateMismatchAndClipRate(reads: Iterable[AlignmentRecord],
                                                      reference: String,
                                                      referenceStartPos: Long): (Double, Double, Double) = {
    val (bases, mismatchedBases, clippedBases) = reads.map(r => {
      val cigar: List[CigarElement] = CIGAR_CODEC.decode(r.getCigar).getCigarElements
      val readSequence = r.getSequence
      var refIdx = (r.getStart - referenceStartPos).toInt
      var readIdx = 0
      var mismatches = 0
      var clips = 0

      // loop over cigar elements
      cigar.foreach(i => {
        i.getOperator match {
          case CigarOperator.M | CigarOperator.X | CigarOperator.EQ => {
            (0 until (i.getLength - 1)).foreach(j => {
              // do we have a mismatch
              if (readSequence(readIdx) != reference(refIdx)) {
                mismatches += 1
              }
              readIdx += 1
              refIdx += 1
            })
          }
          case CigarOperator.S => {
            // we've got a clip
            clips += i.getLength
            readIdx += i.getLength
          }
          case _ => {
            if (i.getOperator.consumesReferenceBases) {
              refIdx += i.getLength
            }
            if (i.getOperator.consumesReadBases) {
              readIdx += i.getLength
            }
          }
        }
      })

      (readIdx + 1, mismatches, clips)
    }).reduce((t1: (Int, Int, Int), t2: (Int, Int, Int)) => {
      (t1._1 + t2._1, t1._2 + t2._2, t1._3 + t1._3)
    })

    (mismatchedBases.toDouble / bases.toDouble, clippedBases.toDouble / bases.toDouble, bases.toDouble / reference.length.toDouble)
  }
}

import ReassemblyExplorer._

class ReassemblyExplorer(kmerLength: Int,
                         reference: RDD[NucleotideContigFragment],
                         sd: SequenceDictionary,
                         contigLengths: Map[String, Long],
                         activeRegionHighCoverageThreshold: Double,
                         activeRegionLowCoverageThreshold: Double,
                         activeRegionMismatchRateThreshold: Double,
                         activeRegionClipRateThreshold: Double) extends Explorer with Logging {

  val totalAssembledReferenceLength = contigLengths.values.sum

  val companion: ExplorerCompanion = ReassemblyExplorer

  def discoverRegion(rr: (Iterable[AlignmentRecord], NucleotideContigFragment)): Iterable[Observation] = RegionDiscovery.time {
    val (reads, reference) = rr

    // extract the coordinates and sequence from the reference
    val coordinates = ReferenceRegion(reference).get
    val refSeq = reference.getFragmentSequence.toUpperCase

    // compute activity statistics
    val (mismatchRate, clipRate, coverage) = CheckActivity.time {
      ReassemblyExplorer.calculateMismatchAndClipRate(reads,
        refSeq,
        coordinates.start)
    }

    // is the region that we're currently looking at active?
    if ((mismatchRate > activeRegionMismatchRateThreshold ||
      clipRate > activeRegionClipRateThreshold) &&
      (coverage > activeRegionLowCoverageThreshold &&
        coverage < activeRegionHighCoverageThreshold)) {
      ReassemblingRegion.time {
        log.info("Reassembling active region " + coordinates)

        // reassemble us up a region
        val graphs = BuildingGraph.time {
          KmerGraph(kmerLength,
            Seq((coordinates, refSeq)),
            reads.toSeq)
        }

        // turn the reassembly graph into observations
        ObservingGraph.time {
          graphs.flatMap(_.toObservations)
        }
      }
    } else {
      InactiveReads.time {
        log.info("Observing inactive region " + coordinates)
        var refPos = coordinates.start
        var readId = refPos.toLong << 32
        reads.flatMap(r => {
          readId += 1L
          ReadExplorer.readToObservations(r, readId)
        }) ++ refSeq.map(base => {
          val observation = new Observation(ReferencePosition(coordinates.referenceName, refPos),
            base.toString)

          // increment site
          refPos += 1

          observation
        })
      }
    }
  }

  def discover(reads: RDD[AlignmentRecord]): RDD[Observation] = {
    val joinReadsAndContigs = JoiningReads.time {
      // zip reference contig fragments with uuids, cache
      val refIds = reference.zipWithUniqueId()
        .map(vk => (vk._2, vk._1))
        .cache()

      // filter mapped reads, join with reference contigs, then extract contig ids
      val joinWithId = ShuffleRegionJoin.partitionAndJoin(reference.context,
        refIds,
        reads.filter(_.getReadMapped),
        sd,
        totalAssembledReferenceLength / reads.partitions.size)
        .flatMap(kv => {
          val ((fragmentId, fragment), read) = kv

          val fStart = fragment.getFragmentStartPosition
          val fEnd = fStart + fragment.getFragmentSequence.length

          if (fStart <= read.getStart && fEnd >= read.getEnd) {
            Some((kv._1._1, kv._2))
          } else {
            None
          }
        })

      // merge together reads that are from the same contig fragment, then
      // join against reference contigs, and throw away id
      // ideally, we should optimize the partitioning; this is a to-do for later
      val jrdd = joinWithId.groupByKey()
        .join(refIds)
        .map(kv => kv._2)

      // we are done with the original reference contig rdd, so we can unpersist
      refIds.unpersist()

      jrdd
    }

    // reassemble our regions into observations
    ProcessingRegions.time {
      joinReadsAndContigs.flatMap(discoverRegion)
    }
  }
}
