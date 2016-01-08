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
import org.bdgenomics.adam.rich.{ DecadentRead, RichAlignmentRecord }
import org.bdgenomics.adam.util.MdTag
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

  def mdTagToMismatchPositions(mdTag: MdTag, cigar: List[CigarElement]): Seq[Int] = {
    var idx = 0
    val insertions = cigar.map(c => {
      (c, c.getLength)
    }).map(kv => {
      val r = (kv._1, idx)
      idx += kv._2
      r
    }).flatMap(kv => {
      val (ce, i) = kv
      if (ce.getOperator == CigarOperator.I) {
        (0 until ce.getLength).map(_ + i)
      } else {
        Seq.empty
      }
    })

    val deletions = mdTag.deletions
    val oriPositions = mdTag.mismatches.keys
    var mismatchPositions = oriPositions.zip(oriPositions)
    for (iPos <- insertions) {
      mismatchPositions = mismatchPositions.map({ case (p, i) => if (i <= iPos) (p + 1, i) else (p, i) })
    }
    for ((dPos, _) <- deletions) {
      mismatchPositions = mismatchPositions.map({ case (p, i) => if (i > dPos) (p - 1, i) else (p, i) })
    }
    mismatchPositions.map({ case (p, i) => p.toInt }).toSeq
  }

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

    val firstOfPair = read.getReadNum == 0
    // get cigar, md tag, and phred scores for bases
    val cigar: List[CigarElement] = richRead.samtoolsCigar.getCigarElements
    val quals = richRead.qualityScores
    val mdString = read.getMismatchingPositions
    val mismatchPositions: Option[Seq[Int]] = if (mdString != null && mdString != "")
      Some(mdTagToMismatchPositions(MdTag(read.getMismatchingPositions,
        if (cigar.head.getOperator == CigarOperator.S) cigar.head.getLength else 0,
        richRead.samtoolsCigar), cigar))
    else None

    // observations
    var observations = Seq[Observation]()

    // position in the current read
    var readPos = 0

    // get the sum of mismatching bases
    val qscores: Option[Seq[Int]] = mismatchPositions.map(l => {
      l.map(p => {
        quals(p)
      })
    })

    val mismatchQScoreSum = qscores.map(_.sum)

    // Helper function to get the unclipped read length (hard or soft) from CIGAR
    def unclippedLenFromCigar(cigar: Cigar): Int = {
      cigar.getCigarElements.map(ce => ce.getOperator match {
        case CigarOperator.D | CigarOperator.N | CigarOperator.P => 0
        case _ => ce.getLength
      }).sum
    }

    def alignedLenFromCigar(cigar: Cigar): Int = {
      cigar.getCigarElements.map(ce => ce.getOperator match {
        case CigarOperator.D | CigarOperator.N | CigarOperator.P | CigarOperator.H | CigarOperator.S => 0
        case _ => ce.getLength
      }).sum
    }

    // Helper function to calculate the length of an element, if it is a clipping element
    def basesTrimmed(cigarElement: CigarElement): Int = {
      cigarElement.getOperator match {
        case CigarOperator.S | CigarOperator.H => cigarElement.getLength
        case _                                 => 0
      }
    }

    // Set up variables to help with tracking the distance from indels, and
    // the distance from the current allele to the first and last trimmed base
    // within this read.
    val readLen = unclippedLenFromCigar(richRead.samtoolsCigar)
    val trimmedFromStart = basesTrimmed(cigar.head)
    val trimmedFromEnd = basesTrimmed(cigar.last)
    var softclippedBases = 0
    val alignedLen = alignedLenFromCigar(richRead.samtoolsCigar)

    val cigarLenOps = cigar.zipWithIndex.map({
      case (ce: CigarElement, idx: Int) =>
        (idx, (ce.getLength, ce.getOperator))
    }).toMap

    val insertions = cigarLenOps.filter({ case (idx, (len, op)) => op == CigarOperator.I })
    val deletions = cigarLenOps.filter({ case (idx, (len, op)) => op == CigarOperator.D })

    def makeNaieveDistanceVec(idx: Int, len: Int, del: Boolean): Vector[Int] = {
      val lpre = (0 until idx).map(cigarLenOps(_)._1).sum
      val lpost = ((idx + 1) until cigar.length).map(cigarLenOps(_)._1).sum
      ((1 to lpre).reverse ++ (if (del) Vector.empty[Int] else Vector.fill(len)(0)) ++ (1 to lpost).map(-_)).toVector
    }

    val insertionDistVecs = insertions.map({ case (idx, (len, _)) => makeNaieveDistanceVec(idx, len, false) })
    val deletionDistVecs = deletions.map({ case (idx, (len, _)) => makeNaieveDistanceVec(idx, len, true) })

    val posToInsDist: Option[Vector[Int]] = if (insertionDistVecs.size > 0) Some(insertionDistVecs.transpose.map(l => l.minBy(Math.abs(_))).toVector) else None
    val posToDelDist: Option[Vector[Int]] = if (deletionDistVecs.size > 0) Some(deletionDistVecs.transpose.map(l => l.minBy(Math.abs(_))).toVector) else None

    def getTags(read: RichAlignmentRecord): Option[Seq[org.bdgenomics.adam.models.Attribute]] = {
      try {
        Option(read.tags)
      } catch {
        case e: NullPointerException => None
      }
    }

    val tags: Option[Seq[org.bdgenomics.adam.models.Attribute]] = getTags(richRead)
    val mateRescue: Boolean = tags.getOrElse(Seq()).exists(a => a.tag == "XT" && a.value == "M")

    def processAlignmentMatch() {
      observations = AlleleObservation(ReferencePosition(contig, pos),
        1,
        sequence(readPos).toString,
        quals(readPos),
        mapq,
        negativeStrand,
        firstOfPair,
        readPos - softclippedBases,
        alignedLen,
        posToInsDist.map(_(readPos)),
        posToDelDist.map(_(readPos)),
        trimmedFromStart,
        trimmedFromEnd,
        readLen,
        mismatchQScoreSum,
        mateRescue,
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
          readPos - softclippedBases,
          alignedLen,
          posToInsDist.map(_(readPos)),
          posToDelDist.map(_(readPos)),
          trimmedFromStart,
          trimmedFromEnd,
          readLen,
          mismatchQScoreSum,
          mateRescue,
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
          readPos - softclippedBases,
          alignedLen,
          posToInsDist.map(_(readPos)),
          posToDelDist.map(_(readPos)),
          trimmedFromStart,
          trimmedFromEnd,
          readLen,
          mismatchQScoreSum,
          mateRescue,
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
          val clippedLen = cigar(i).getLength
          readPos += clippedLen
          softclippedBases += clippedLen
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
