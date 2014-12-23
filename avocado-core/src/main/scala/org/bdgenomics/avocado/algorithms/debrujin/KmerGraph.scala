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

import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.adam.models.{ ReferencePosition, ReferenceRegion }
import org.bdgenomics.adam.rich.RichAlignmentRecord
import org.bdgenomics.adam.util.PhredUtils
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
          newKmer.phred = PhredUtils.phredToSuccessProbability(q) :: newKmer.phred
          newKmer.mapq = PhredUtils.phredToSuccessProbability(mapq) :: newKmer.mapq

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
          var newKmer: Kmer = null
          val phred = List(PhredUtils.phredToSuccessProbability(q))
          val mapQ = List(PhredUtils.phredToSuccessProbability(mapq))

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

    val graph = new KmerGraph(kmerMap.values.toArray)

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
case class KmerGraph(protected val kmers: Array[Kmer]) {

  // source/sink kmers
  private def sourceKmers = kmers.filter(k => k.predecessors.length == 0 && k.refPos.isDefined)
  private def sinkKmers = kmers.filter(k => k.successors.length == 0 && k.refPos.isDefined)

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

  def size: Int = kmers.length

  def nonRefSize: Int = kmers.filter(_.refPos.isEmpty).length

  def sources: Int = sourceKmers.length

  def sinks: Int = sinkKmers.length
}
