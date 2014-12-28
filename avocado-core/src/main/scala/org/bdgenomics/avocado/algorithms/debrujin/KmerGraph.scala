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
package org.bdgenomics.avocado.algorithms.debrujin

import org.bdgenomics.adam.models.{ ReferencePosition, ReferenceRegion }
import org.bdgenomics.adam.rich.RichAlignmentRecord
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.formats.avro.AlignmentRecord
import scala.annotation.tailrec
import scala.collection.mutable.{ HashMap }

object KmerGraph {

  /**
   * Creates a new de Brujin Graph that shows connections between k-mers. k-mers are inserted
   * from a read dataset that is provided.
   *
   * @param kmerLength Kmer length.
   * @param readLength Read length.
   * @param regionLength Length of the active region.
   * @param reference Reference sequence.
   * @param reads Reads to insert into graph.
   * @param flankLength Length of flanking sequence to take from reference sequence.
   * @param maxEntries Maximum number of times a loop can be entered.
   * @param removeSpurs If true, removes spurs (low quality paths) from the graph.
   * @param lowCoverageTrimmingThreshold Threshold for trimming nodes that are not covered by many reads.
   * @return Returns a new de Brujin graph.
   */
  def apply(kmerLength: Int,
            references: Seq[(ReferenceRegion, String)],
            reads: Seq[RichAlignmentRecord],
            removeSpurs: Boolean = false): KmerGraph = {

    assert(references.length == 1, "For now, only support a single reference in the graph.")

    val kmerMap = HashMap[String, Kmer]()

    @tailrec def addReferenceKmers(iter: Iterator[String],
                                   pos: ReferencePosition,
                                   lastKmer: Kmer) {
      // do we have k-mers left?
      if (iter.hasNext) {
        var newKmer: Kmer = null
        val ks = iter.next

        // if we have a predecessor, populate the predecessor fields
        if (lastKmer != null) {
          val kl = List(lastKmer)
          newKmer = Kmer(ks, Some(pos), predecessors = kl)
          lastKmer.successors = newKmer :: lastKmer.successors
        } else {
          newKmer = Kmer(ks, Some(pos))
        }

        // add the kmer to the graph
        kmerMap(ks) = newKmer

        // update the position for the next k-mer
        val newPos = ReferencePosition(pos.referenceName, pos.pos + 1L)

        addReferenceKmers(iter, newPos, newKmer)
      }
    }

    references.foreach(p => {
      val (region, reference) = p
      addReferenceKmers(reference.sliding(kmerLength),
        ReferencePosition(region.referenceName, region.start),
        null)
    })

    @tailrec def addReadKmers(kmerIter: Iterator[String],
                              qualIter: Iterator[Char],
                              mapq: Int,
                              lastKmer: Kmer) {
      // do we have k-mers left?
      if (kmerIter.hasNext) {
        assert(qualIter.hasNext, "Should still have qualities as long as we've got k-mers.")
        var newKmer: Kmer = null
        val ks = kmerIter.next
        val q = qualIter.next

        // have we already seen this k-mer?
        if (kmerMap.contains(ks)) {
          newKmer = kmerMap(ks)

          // add phred score and mapq
          newKmer.phred = q.toInt :: newKmer.phred
          newKmer.mapq = mapq :: newKmer.mapq

          // do we have a predecessor? if so, perform book keeping...
          if (lastKmer != null) {
            // if we have a predecessor, and it isn't in the k-mer map, then add it
            if (newKmer.predecessors.filter(_.kmerSeq == lastKmer.kmerSeq).length == 0) {
              newKmer.predecessors = lastKmer :: newKmer.predecessors
            }

            // if this k-mer isn't in the successor list, then add it
            if (lastKmer.successors.filter(_.kmerSeq == ks).length == 0) {
              lastKmer.successors = newKmer :: lastKmer.successors
            }
          }
        } else {
          val phred = List(q.toInt)
          val mapQ = List(mapq)

          // if we have a predecessor, populate the predecessor fields
          if (lastKmer != null) {
            val kl = List(lastKmer)
            newKmer = Kmer(ks, phred = phred, mapq = mapQ, predecessors = kl)
            lastKmer.successors = newKmer :: lastKmer.successors
          } else {
            newKmer = Kmer(ks, phred = phred, mapq = mapQ)
          }

          // add the kmer to the graph
          kmerMap(ks) = newKmer
        }

        addReadKmers(kmerIter, qualIter, mapq, newKmer)
      }
    }

    // loop over reads and collect statistics
    reads.foreach(r => {
      addReadKmers(r.getSequence.toString.sliding(kmerLength),
        r.getQual.toString.toIterator,
        r.getMapq,
        null)
    })

    val graph = new KmerGraph(kmerMap.values.toArray, kmerLength)

    if (removeSpurs) {
      graph.removeSpurs()
    }

    graph
  }
}

/**
 * Graph showing connections between kmers.
 *
 * @param kmers k-mers which were observed.
 */
case class KmerGraph(protected val kmers: Array[Kmer],
                     protected val kmerLength: Int) {

  // source/sink kmers
  private val allSourceKmers = kmers.filter(_.predecessors.length == 0)
  private val sourceKmers = allSourceKmers.filter(_.refPos.isDefined)
  private val allSinkKmers = kmers.filter(_.successors.length == 0)
  private val sinkKmers = allSinkKmers.filter(_.refPos.isDefined)

  /**
   * Removes spurs from graph. Spurs are segments of the graph that do not connect to
   * the source or sink of the graph. While these spurs do not contribute spurious haplotypes,
   * they make the haplotype enumeration process more expensive.
   */
  def removeSpurs() {
    ???
  }

  override def toString(): String = {
    "Sources: " + sourceKmers.map(_.kmerSeq).reduce(_ + ", " + _) + "\n" +
      "Sinks: " + sinkKmers.map(_.kmerSeq).reduce(_ + ", " + _) + "\n" +
      kmers.map(_.toString).reduce(_ + "\n" + _)
  }

  /**
   * Prints this de Brujin graph out in Graphviz format. This can be used with the Dot software
   * to visualize graph creation.
   *
   * @return Returns a string describing this de Brujin graph as a directed Graphviz graph.
   */
  def toDot(): String = {
    "digraph kg { \n" +
      kmers.map(_.toDot).reduce(_ + "\n" + _) + "\n}\n"
  }

  /**
   * Converts this graph of k-mers into a set of observations. The observations
   * include both the reference threading as well as the actual read derived observations
   * that can be used for genotyping. This is performed without realignment.
   *
   * @return Returns a seq of observations.
   */
  def toObservations: Seq[Observation] = {
    def buildReadObservations(kmer: Kmer,
                              pos: ReferencePosition,
                              length: Int,
                              allele: String): Seq[Observation] = {
      (0 until kmer.multiplicity).map(i => {
        AlleleObservation(pos,
          length,
          allele,
          kmer.phred(i),
          kmer.mapq(i),
          false,
          "")
      })
    }

    def buildReferenceObservations(kmer: Kmer): Seq[Observation] = {
      val site = Seq(new Observation(kmer.refPos.get, kmer.kmerSeq.take(1)))
      // did we have any reads covering this site?
      // if not, just emit that we observed the site
      if (kmer.multiplicity == 0) {
        site
      } else {
        site ++ buildReadObservations(kmer, kmer.refPos.get, 1, kmer.kmerSeq.take(1))
      }
    }

    @tailrec def crawl(nextKmer: Kmer,
                       observations: Seq[Observation],
                       pending: List[Kmer],
                       allele: String,
                       branchPoint: ReferencePosition,
                       branches: List[(ReferencePosition, Kmer)]): Seq[Observation] = {

      if (nextKmer == null) {
        observations
      } else {
        val (newNextKmer,
          newObservations,
          nowPending,
          newAllele,
          newBranch,
          newBranches): (Kmer, Seq[Observation], List[Kmer], String, ReferencePosition, List[(ReferencePosition, Kmer)]) = if (branchPoint != null && (nextKmer.refPos.isEmpty || !pending.isEmpty)) {
          // if the current k-mer is unmapped, then we're building an alt
          // else, we're cleaning up an alt
          if (nextKmer.refPos.isEmpty) {
            // extend allele and pending list
            val newAllele = if (allele != null) {
              allele + nextKmer.kmerSeq.take(1)
            } else {
              nextKmer.kmerSeq.take(1)
            }
            val newPending = nextKmer :: pending

            // for now, we require alts to not diverge
            assert(nextKmer.successors.length == 1,
              "Unexpected divergence at: " + nextKmer.toDetailedString)
            (nextKmer.successors.head,
              observations,
              newPending,
              newAllele,
              branchPoint,
              branches)
          } else {
            // where have we connected back to?
            val refSink = nextKmer.refPos.get

            // reverse pending k-mers
            val reversed = pending.reverse

            // create observations of the alt allele
            val altObs = reversed.drop(kmerLength - 1)
              .map(buildReadObservations(_,
                ReferencePosition(branchPoint.referenceName,
                  branchPoint.pos + kmerLength),
                (refSink.pos - branchPoint.pos).toInt - kmerLength + 1,
                allele.drop(kmerLength - 1)))
              .reduce(_ ++ _)

            // now, count (k - 1) from the start and build reference observations
            var refPoint = ReferencePosition(branchPoint.referenceName,
              branchPoint.pos + 1)
            val refObs = reversed.take(kmerLength - 1)
              .map(k => {
                // take observations
                val obs = buildReadObservations(k,
                  refPoint,
                  1,
                  k.kmerSeq.take(1))

                // increment position
                refPoint = ReferencePosition(refPoint.referenceName,
                  refPoint.pos + 1)

                // emit observations
                obs
              }).reduce(_ ++ _)

            // combine observations into a new observation set and recurse
            val newObs = observations ++ altObs ++ refObs
            (nextKmer,
              newObs,
              List(),
              null,
              null,
              branches)
          }
        } else {
          assert(nextKmer.refPos.isDefined,
            "k-mer must be mapped.")
          val pos = nextKmer.refPos.get
          assert(branchPoint == null ||
            pos.referenceName == branchPoint.referenceName &&
            pos.pos == branchPoint.pos + 1,
            "Two reference k-mers must be adjacent in reference space.")

          // build observations
          val newObservations = buildReferenceObservations(nextKmer)

          // if we have a successor, we'll take that and continue on
          // else, we'll move to the next branch
          val (next, newBranches) = if (nextKmer.successors.isEmpty && branches.isEmpty) {
            ((null, null), branches)
          } else if (nextKmer.successors.isEmpty) {
            (branches.head, branches.drop(1))
          } else {
            ((null, nextKmer.successors.filter(_.refPos.isDefined).head),
              nextKmer.successors.filter(_.refPos.isEmpty).map(v => (pos, v)) ::: branches)
          }

          (next._2,
            observations ++ newObservations,
            pending,
            null,
            next._1,
            newBranches)
        }

        // recurse and compute next iteration
        crawl(newNextKmer,
          newObservations,
          nowPending,
          newAllele,
          newBranch,
          newBranches)
      }
    }

    crawl(sourceKmers.head,
      Seq(),
      List(),
      null,
      null,
      sourceKmers.drop(1).map(p => (null.asInstanceOf[ReferencePosition], p)).toList)
  }

  def size: Int = kmers.length

  def nonRefSize: Int = kmers.filter(_.refPos.isEmpty).length

  def sources: Int = sourceKmers.length

  def sinks: Int = sinkKmers.length
}
