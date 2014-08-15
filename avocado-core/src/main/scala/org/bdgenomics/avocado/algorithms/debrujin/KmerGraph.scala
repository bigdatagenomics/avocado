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
import org.bdgenomics.adam.rich.RichAlignmentRecord
import scala.annotation.tailrec
import scala.collection.mutable.{ HashSet, HashMap, SortedSet }

object KmerGraph {

  /**
   * Chops a read into k-mers.
   *
   * @param read Read to chop into k-mers.
   * @param k _k_-mer length
   * @return Returns a seq of strings with length _k_.
   */
  private def getKmerSequences(read: AlignmentRecord, k: Int): Seq[String] = {
    val readSeq: String = read.getSequence.toString
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

  @tailrec private[debrujin] final def mergeGraphs(graph1: Map[String, Set[Kmer]],
                                                   graph2: Map[String, Set[Kmer]]): Map[String, Set[Kmer]] = {
    def mergeKmer(set: Set[Kmer], kmer: Kmer): Set[Kmer] = {
      // get all matching kmers
      val likeKmers = set.filter(_.kmerSeq == kmer.kmerSeq)

      if (likeKmers.size >= 1) {
        // are any of the kmers in flanking sequence?
        val inFlank = likeKmers.map(_.inFlank).reduce(_ || _) || kmer.inFlank

        // update multiplicity of kmer
        set.filter(_.kmerSeq != kmer.kmerSeq) + new Kmer(kmer.kmerSeq,
          likeKmers.map(_.weight).sum,
          inFlank)
      } else {
        set + kmer
      }
    }

    @tailrec def mergeKmerSet(prefix: String,
                              set1: Set[Kmer],
                              set2: Set[Kmer]): (String, Set[Kmer]) = {
      if (set2.isEmpty) {
        (prefix, set1)
      } else {
        val newSet1 = mergeKmer(set1, set2.head)
        mergeKmerSet(prefix, newSet1, set2.drop(1))
      }
    }

    def mergePrefix(graph: Map[String, Set[Kmer]],
                    node: (String, Set[Kmer])): Map[String, Set[Kmer]] = {
      val (prefix, set) = node

      // if prefix in graph, merge, else add
      if (graph.contains(prefix)) {
        graph.filterKeys(_ != node) + mergeKmerSet(prefix, graph(prefix), set)
      } else {
        graph + node
      }
    }

    if (graph2.isEmpty) {
      graph1
    } else {
      val merged = mergePrefix(graph1, graph2.head)
      mergeGraphs(merged, graph2.drop(1))
    }
  }

  private[debrujin] def buildPrefixMap(kmerSequences: Seq[String],
                                       isFlank: Boolean = false): Map[String, Set[Kmer]] = {
    val kmers = kmerSequences.groupBy(x => x) // group by sequence
      .map(kv => Kmer(kv._1, kv._2.size, isFlank)) // generate kmer w/ sequence and count
      .groupBy(kmer => kmer.prefix)

    kmers.map(kv => (kv._1, kv._2.toSet))
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
   * @param maxEntries Maximum number of times a loop can be entered.
   * @param removeSpurs If true, removes spurs (low quality paths) from the graph.
   * @param lowCoverageTrimmingThreshold Threshold for trimming nodes that are not covered by many reads.
   * @return Returns a new de Brujin graph.
   */
  def apply(kmerLength: Int,
            readLength: Int,
            regionLength: Int,
            reference: String,
            reads: Seq[RichAlignmentRecord],
            flankLength: Int,
            maxEntries: Int = 5,
            removeSpurs: Boolean = false,
            lowCoverageTrimmingThreshold: Option[Double] = None): KmerGraph = {

    val graph = new KmerGraph(kmerLength, readLength, regionLength, reference, flankLength, maxEntries)
    graph.insertReads(reads)
    if (removeSpurs || lowCoverageTrimmingThreshold.isDefined) {
      lowCoverageTrimmingThreshold.foreach(t => graph.trimLowCoverageKmers(t))
      graph.removeSpurs()
    }
    graph
  }

  def apply(kmerLength: Int,
            readLength: Int,
            regionLength: Int,
            reference: String,
            flankLength: Int,
            maxEntries: Int): KmerGraph = {
    new KmerGraph(kmerLength, readLength, regionLength, reference, flankLength, maxEntries)
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
class KmerGraph(val kmerLength: Int,
                val readLength: Int,
                val regionLength: Int,
                val reference: String,
                flankLength: Int,
                maxEntries: Int) {

  assert(flankLength >= kmerLength, "Flanking sequence length must be longer than the k-mer length.")

  // source/sink kmerLength
  val sourceKmer = new Kmer(reference.take(kmerLength))
  val sinkKmer = new Kmer(reference.takeRight(kmerLength))

  // The actual kmer graph consists of unique K-1 prefixes and kmers connected
  // by vertices.
  private var kmerGraph = Map[String, Set[Kmer]]()

  // add flanking sequences
  addFlanks(reference, flankLength)

  def prefixSet = kmerGraph.keys.toSet

  def kmers = kmerGraph.flatMap(_._2)

  def kmerSequences = kmers.map(_.kmerSeq)

  // Paths through the kmer graph are in order of decreasing total mult.
  val maxPathDepth = regionLength + readLength - kmerLength + 1
  lazy val allPaths = enumerateAllPaths(maxPathDepth, maxEntries)

  def insertReads(reads: Seq[RichAlignmentRecord]) = {
    val kmerSequences = reads.flatMap(read => KmerGraph.getKmerSequences(read, kmerLength))
    addSequencesToGraph(kmerSequences)
  }

  def insertRead(read: AlignmentRecord) = {
    val kmerSequences = KmerGraph.getKmerSequences(read, kmerLength)
    addSequencesToGraph(kmerSequences)
  }

  def addFlanks(reference: String, flankLength: Int) {
    val startFlank = reference.take(flankLength)
    val endFlank = reference.takeRight(flankLength)
    val kmerSequences = KmerGraph.getKmerSequences(startFlank, kmerLength) ++ KmerGraph.getKmerSequences(endFlank, kmerLength)
    addSequencesToGraph(kmerSequences, true)
  }

  def addSequencesToGraph(kmerSequences: Seq[String], isFlank: Boolean = false) {
    kmerGraph = KmerGraph.mergeGraphs(kmerGraph, KmerGraph.buildPrefixMap(kmerSequences, isFlank))
  }

  /**
   * Performs depth first search and scores all paths through the final graph.
   */
  def enumerateAllPaths(maxDepth: Int, maxEntries: Int): SortedSet[KmerPath] = {
    val paths = SortedSet[KmerPath]()
    val hitCount = HashMap[String, Int]()

    def pathToSink(node: Kmer, currentPath: Seq[Kmer], depth: Int): Unit = {
      val seq = node.kmerSeq
      // how many times have we gone into this kmer?
      val hits = hitCount.getOrElse(seq, 0)

      if (node.nextPrefix.equals(sinkKmer.nextPrefix)) {
        // if we've hit the end, add our path to the list
        val path = new KmerPath(currentPath :+ node)
        paths.add(path)
      } else if (depth <= maxDepth && hits <= maxEntries) {
        // entering kmer
        hitCount(seq) = hits + 1

        // walk next kmers
        val nextNodes = kmerGraph.get(node.nextPrefix)

        if (nextNodes.isDefined && nextNodes.nonEmpty) {
          nextNodes.get.map(next => pathToSink(next, currentPath :+ node, depth + 1))
        }

        // exiting kmer
        hitCount(seq) = hits
      }
    }

    pathToSink(sourceKmer, Seq.empty, 0)
    paths
  }

  /**
   * Removes low coverage kmers from the graph. This is used to remove graph bubbles that do
   * not have significant evidence, which limits the number of haplotypes created. After this
   * is called, spur removal should be run. Removes kmers that have coverage lower than a given
   * ratio times the median kmer coverage. Median is chosen instead of mean to avoid biasing
   * versus repeats.
   *
   * @param ratio Ratio used to set coverage threshold.
   *
   * @see removeSpurs
   */
  def trimLowCoverageKmers(ratio: Double) {
    // get median coverage of kmers that are covered by more than 1 read
    val kmerWeights = kmers.map(_.weight).filter(_ > 1).toSeq.sortWith(_ < _)
    val median = kmerWeights(kmerWeights.length / 2)

    // threshold for trimming
    val threshold = (median.toDouble * ratio).toInt

    kmerGraph = kmerGraph.mapValues(s => s.filter(k => k.weight >= threshold || k.inFlank))
      .filter(kv => !kv._2.isEmpty)
  }

  /**
   * Removes spurs from graph. Spurs are segments of the graph that do not connect to
   * the source or sink of the graph. While these spurs do not contribute spurious haplotypes,
   * they make the haplotype enumeration process more expensive.
   */
  @tailrec final def removeSpurs() {
    def count: Int = kmerGraph.map(kv => kv._2.size).sum

    val lastCount = count

    val searchPrefixes = kmers.map(_.nextPrefix).toSeq
    val sourceSearchPrefix = sourceKmer.nextPrefix
    val sourcePrefix = sourceKmer.prefix
    val sinkPrefix = sinkKmer.prefix
    val sinkSearchPrefix = sinkKmer.nextPrefix
    val prefixes = prefixSet

    kmerGraph = kmerGraph.filter(kv => {
      val (prefix, nodeKmers) = kv

      // filter out kmers who do not have a predecessor
      searchPrefixes.contains(prefix) || prefix == sourceSearchPrefix ||
        prefix == sourcePrefix || prefix == sinkPrefix
    }).mapValues(ks => {
      ks.filter(k => {
        val next = k.nextPrefix

        // filter out kmers who do not point at any other kmers
        prefixes.contains(next) || next == sinkPrefix || next == sinkSearchPrefix
      })
    }).filter(kv => !kv._2.isEmpty)

    // if we removed spurs this iteration, loop and remove spurs
    if (count != lastCount) {
      removeSpurs()
    }
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
