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

import htsjdk.samtools.{ Cigar, CigarElement, CigarOperator }
import org.apache.commons.configuration.{ HierarchicalConfiguration, SubnodeConfiguration }
import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferencePosition
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rich.RichAlignmentRecord
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.AlignmentRecord

object ReadExplorer extends ExplorerCompanion {

  val explorerName: String = "ReadExplorer"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Explorer = {
    new ReadExplorer(stats.referenceObservations)
  }
}

class ReadExplorer(referenceObservations: RDD[Observation]) extends Explorer with Logging {

  val companion: ExplorerCompanion = ReadExplorer

  def readToObservations(r: (AlignmentRecord, Long)): Seq[Observation] = ExploringRead.time {
    val (read, readId) = r
    val richRead: RichAlignmentRecord = RichAlignmentRecord(read)

    // get read start, contig, strand, sample, mapq, and sequence
    var pos: Long = read.getStart
    val contig: String = read.getContig.getContigName
    val negativeStrand: Boolean = read.getReadNegativeStrand
    val sample: String = read.getRecordGroupSample.toString
    // you can't use Option(read.getMapq), because scala calls scala.Predef$.Integer2int...
    val mapq: Option[Int] = if (read.getMapq == null) {
      None
    } else {
      Some(read.getMapq)
    }
    val sequence: String = read.getSequence

    val firstOfPair = read.getFirstOfPair
    // get cigar, md tag, and phred scores for bases
    val cigar: List[CigarElement] = richRead.samtoolsCigar.getCigarElements
    val quals = richRead.qualityScores
    val mdTag = richRead.mdTag

    // observations
    var observations = Seq[Observation]()

    // position in the current read
    var readPos = 0

    def processAlignmentMatch() {
      observations = AlleleObservation(ReferencePosition(contig, pos),
        1,
        sequence(readPos).toString,
        quals(readPos),
        mapq,
        negativeStrand,
        firstOfPair,
        readPos,
        sample,
        readId) +: observations
      readPos += 1
      pos += 1
    }

    def processAlignmentMatchWithLookahead(idx: Int) {
      if (idx + 1 < cigar.length && cigar(idx + 1).getOperator == CigarOperator.I) {
        // the allele includes the matching base
        val alleleLength = 1 + cigar(idx + 1).getLength

        // add the observation
        observations = AlleleObservation(ReferencePosition(contig, pos),
          1,
          sequence.drop(readPos).take(alleleLength),
          quals.drop(readPos).take(alleleLength).reduce(_ + _) / alleleLength,
          mapq,
          negativeStrand,
          firstOfPair,
          readPos,
          sample,
          readId).asInstanceOf[Observation] +: observations

        // increment read pointers
        readPos += alleleLength
        pos += 1
      } else if (idx + 1 < cigar.length && cigar(idx + 1).getOperator == CigarOperator.D) {
        // the allele includes the matching base
        val alleleLength = 1 + cigar(idx + 1).getLength

        // add the observation
        observations = AlleleObservation(ReferencePosition(contig, pos),
          alleleLength,
          sequence.drop(readPos).take(1),
          quals.drop(readPos).head,
          mapq,
          negativeStrand,
          firstOfPair,
          readPos,
          sample,
          readId).asInstanceOf[Observation] +: observations

        // increment read pointers
        readPos += 1
        pos += alleleLength
      } else {
        processAlignmentMatch()
      }
    }

    // loop over cigar elements
    (0 until cigar.length).foreach(i => {
      cigar(i).getOperator match {
        case CigarOperator.M | CigarOperator.X | CigarOperator.EQ => {
          // loop over cigar, but stop before the last base
          // we must look ahead for an insert at the last base
          (0 until (cigar(i).getLength - 1)).foreach(j => {
            // process match
            processAlignmentMatch()
          })

          // look ahead for indel and process
          processAlignmentMatchWithLookahead(i)
        }
        case CigarOperator.I => {
          // no op; handle inserts by looking ahead from match/mismatch operator
        }
        case CigarOperator.D => {
          // no op; handle inserts by looking ahead from match/mismatch operator
        }
        case CigarOperator.S => {
          readPos += cigar(i).getLength
        }
        case CigarOperator.H =>
        case _ => {
          if (cigar(i).getOperator.consumesReferenceBases) {
            pos += cigar(i).getLength
          }
          if (cigar(i).getOperator.consumesReadBases) {
            readPos += cigar(i).getLength
            log.warn("Unexpected cigar operator " + cigar(i) + " in: " + read)
          }
        }
      }
    })

    observations
  }

  def discover(reads: RDD[AlignmentRecord]): RDD[Observation] = {
    ExploringReads.time {
      reads.filter(_.getReadMapped)
        .zipWithUniqueId()
        .flatMap(readToObservations)
    } ++ referenceObservations
  }
}
