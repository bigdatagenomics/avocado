/*
 * Copyright (c) 2013-2014. Regents of the University of California
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

package org.bdgenomics.avocado.algorithms.debrujin

import org.bdgenomics.adam.avro.ADAMRecord
import org.bdgenomics.adam.rich.RichADAMRecord
import scala.collection.mutable.{ HashSet }

object KmerGraph {

  private def getKmerSequences(read: ADAMRecord, k: Int): Seq[String] = {
    val readSeq = read.getSequence.toString
    val readLen = readSeq.size
    val offsets = (0 until readLen - k + 1).toSeq
    offsets.map(idx => readSeq.substring(idx, idx + k))
  }

  private def buildPrefixMap(kmerSequences: Seq[String]): Map[String, HashSet[Kmer]] = {
    val kmers = kmerSequences.groupBy(x => x) // group by sequence
      .map(kv => Kmer(kv._1, kv._2.size)) // generate kmer w/ sequence and count
      .groupBy(kmer => kmer.prefix)

    kmers.map(kv => (kv._1, HashSet[Kmer](kv._2.toSeq: _*)))
  }

  def apply(kmerLength: Int, readLen: Int, regionLen: Int, reference: String, reads: Seq[RichADAMRecord], removeSpurs: Boolean = false): KmerGraph = {

    val graph = new KmerGraph(kmerLength, readLen, regionLen, reference)
    graph.insertReads(reads)
    if (removeSpurs) graph.removeSpurs()
    graph
  }

  def apply(kmerLength: Int, readLen: Int, regionLen: Int, reference: String): KmerGraph = {
    new KmerGraph(kmerLength, readLen, regionLen, reference)
  }
}

/**
 * Graph showing connections between kmers.
 *
 * @param kmerLength Kmer length.
 * @param readLen Read length.
 * @param regionLen Length of the active region.
 */
class KmerGraph(val kmerLength: Int, val readLen: Int, val regionLen: Int, val reference: String) {

  // source/sink kmerLength
  val sourceKmer = new Kmer(reference.take(kmerLength))
  val sinkKmer = new Kmer(reference.takeRight(kmerLength))

  // The actual kmer graph consists of unique K-1 prefixes and kmers connected
  // by vertices.
  private val kmerGraph = new scala.collection.mutable.HashMap[String, scala.collection.mutable.HashSet[Kmer]]

  def prefixSet = kmerGraph.keys.toSet

  def kmers = kmerGraph.flatMap(_._2)

  def kmerSequences = kmers.map(_.kmerSeq)

  // Paths through the kmer graph are in order of decreasing total mult.
  val maxPathDepth = regionLen + readLen - kmerLength + 1
  lazy val allPaths = enumerateAllPaths(maxPathDepth)

  def insertReads(reads: Seq[RichADAMRecord]) = {
    val kmerSequences = reads.flatMap(read => KmerGraph.getKmerSequences(read, kmerLength))
    kmerGraph ++= KmerGraph.buildPrefixMap(kmerSequences)
  }

  def insertRead(read: ADAMRecord) = {
    val kmerSequences = KmerGraph.getKmerSequences(read, kmerLength)
    kmerGraph ++= KmerGraph.buildPrefixMap(kmerSequences)
  }

  /**
   * Performs depth first search and scores all paths through the final graph.
   */
  def enumerateAllPaths(maxDepth: Int): collection.mutable.PriorityQueue[KmerPath] = {
    val queue = new collection.mutable.PriorityQueue[KmerPath]()

    def pathToSink(node: Kmer, currentPath: Seq[Kmer], depth: Int): Unit = {
      if (node.nextPrefix.equals(sinkKmer.nextPrefix) || depth > maxDepth) {
        queue.enqueue(new KmerPath(currentPath :+ node))
      }
      val nextNodes = kmerGraph.get(node.nextPrefix)
      if (depth <= maxDepth && nextNodes.isDefined && nextNodes.nonEmpty) {
        nextNodes.get.map(next => pathToSink(next, currentPath :+ node, depth + 1))
      }
    }

    pathToSink(sourceKmer, Seq.empty, 0)
    queue
  }

  def exciseKmer(kmer: Kmer) = {
    val prefixSet = kmerGraph.get(kmer.prefix)
    prefixSet.foreach(set => kmerGraph.put(kmer.prefix, set - kmer))
  }

  /**
   * Removes spurs from graph.
   */
  def removeSpurs(spurThreshold: Int = kmerLength) = {
    // Remove all kmer paths with length <= spurThreshold, connected to the
    // graph source and sink.
    val edgesInSpurCandidate: Set[Kmer] = allPaths.filter(_.len < spurThreshold).flatMap(_.edges).toSet
    val edgePrefixes = edgesInSpurCandidate.groupBy(kmer => kmer.prefix)

    def removeCandidateSpur(prefixWithEdges: (String, Set[Kmer])) = {
      val (prefix, edges) = prefixWithEdges
      val neighbors = kmerGraph(prefix)
      if (neighbors.size - edges.size == 0) {
        kmerGraph.remove(prefix)
      } else {
        edges.foreach(exciseKmer)
      }
    }

    edgePrefixes.foreach(removeCandidateSpur)
  }

}
