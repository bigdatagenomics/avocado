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

  /**
   * Chops a read into k-mers.
   *
   * @param read Read to chop into k-mers.
   * @param k _k_-mer length
   * @return Returns a seq of strings with length _k_.
   */
  private def getKmerSequences(read: ADAMRecord, k: Int): Seq[String] = {
    val readSeq = read.getSequence.toString
    getKmerSequences(readSeq, k)
  }

  /**
   * Chops a generic string into k-mers.
   *
   * @param read Read to chop into k-mers.
   * @param k _k_-mer length
   * @return Returns a seq of strings with length _k_.
   */
  private def getKmerSequences(sequence: String, k: Int): Seq[String] = {
    val seqLen = sequence.size
    val offsets = (0 until seqLen - k + 1).toSeq
    offsets.map(idx => sequence.substring(idx, idx + k))
  }

  private def buildPrefixMap(kmerSequences: Seq[String]): Map[String, HashSet[Kmer]] = {
    val kmers = kmerSequences.groupBy(x => x) // group by sequence
      .map(kv => Kmer(kv._1, kv._2.size)) // generate kmer w/ sequence and count
      .groupBy(kmer => kmer.prefix)

    kmers.map(kv => (kv._1, HashSet[Kmer](kv._2.toSeq: _*)))
  }

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
   * @param removeSpurs If true, removes spurs (low quality paths) from the graph.
   * @return Returns a new de Brujin graph.
   */
  def apply(kmerLength: Int,
            readLength: Int,
            regionLength: Int,
            reference: String,
            reads: Seq[RichADAMRecord],
            flankLength: Int,
            removeSpurs: Boolean = false): KmerGraph = {

    val graph = new KmerGraph(kmerLength, readLength, regionLength, reference, flankLength)
    graph.insertReads(reads)
    if (removeSpurs) graph.removeSpurs()
    graph
  }

  def apply(kmerLength: Int, readLength: Int, regionLength: Int, reference: String, flankLength: Int): KmerGraph = {
    new KmerGraph(kmerLength, readLength, regionLength, reference, flankLength)
  }
}

/**
 * Graph showing connections between kmers.
 *
 * @param kmerLength Kmer length.
 * @param readLength Read length.
 * @param regionLength Length of the active region.
 * @param reference Reference sequence.
 * @param flankLength Length of flanking sequence to take from reference sequence.
 */
class KmerGraph(val kmerLength: Int, val readLength: Int, val regionLength: Int, val reference: String, flankLength: Int) {

  assert(flankLength >= kmerLength, "Flanking sequence length must be longer than the k-mer length.")

  // source/sink kmerLength
  val sourceKmer = new Kmer(reference.take(kmerLength))
  val sinkKmer = new Kmer(reference.takeRight(kmerLength))

  // The actual kmer graph consists of unique K-1 prefixes and kmers connected
  // by vertices.
  private val kmerGraph = new scala.collection.mutable.HashMap[String, scala.collection.mutable.HashSet[Kmer]]

  // add flanking sequences
  addFlanks(reference, flankLength)

  def prefixSet = kmerGraph.keys.toSet

  def kmers = kmerGraph.flatMap(_._2)

  def kmerSequences = kmers.map(_.kmerSeq)

  // Paths through the kmer graph are in order of decreasing total mult.
  val maxPathDepth = regionLength + readLength - kmerLength + 1
  lazy val allPaths = enumerateAllPaths(maxPathDepth)

  def insertReads(reads: Seq[RichADAMRecord]) = {
    val kmerSequences = reads.flatMap(read => KmerGraph.getKmerSequences(read, kmerLength))
    addSequencesToGraph(kmerSequences)
  }

  def insertRead(read: ADAMRecord) = {
    val kmerSequences = KmerGraph.getKmerSequences(read, kmerLength)
    addSequencesToGraph(kmerSequences)
  }

  def addFlanks(reference: String, flankLength: Int) {
    val startFlank = reference.take(flankLength)
    val endFlank = reference.takeRight(flankLength)
    val kmerSequences = KmerGraph.getKmerSequences(startFlank, kmerLength) ++ KmerGraph.getKmerSequences(endFlank, kmerLength)
    addSequencesToGraph(kmerSequences)
  }

  def addSequencesToGraph(kmerSequences: Seq[String]) {
    kmerGraph ++= KmerGraph.buildPrefixMap(kmerSequences)
  }

  /**
   * Performs depth first search and scores all paths through the final graph.
   */
  def enumerateAllPaths(maxDepth: Int): collection.mutable.SortedSet[KmerPath] = {
    val paths = collection.mutable.SortedSet[KmerPath]()

    def pathToSink(node: Kmer, currentPath: Seq[Kmer], depth: Int): Unit = {
      if (node.nextPrefix.equals(sinkKmer.nextPrefix) || depth > maxDepth) {
        val path = new KmerPath(currentPath :+ node)
        paths.add(path)
      }
      val nextNodes = kmerGraph.get(node.nextPrefix)
      if (depth <= maxDepth && nextNodes.isDefined && nextNodes.nonEmpty) {
        nextNodes.get.map(next => pathToSink(next, currentPath :+ node, depth + 1))
      }
    }

    pathToSink(sourceKmer, Seq.empty, 0)
    paths
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

  override def toString(): String = {
    "Source: " + sourceKmer.kmerSeq + "\n" +
      "Sink: " + sinkKmer.kmerSeq + "\n" +
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
      sourceKmer.prefix + " [ shape=box ] ;\n" +
      sinkKmer.nextPrefix + " [ shape=box ] ;\n" +
      kmers.map(_.toDot).reduce(_ + "\n" + _) + "\n}\n"
  }

}
