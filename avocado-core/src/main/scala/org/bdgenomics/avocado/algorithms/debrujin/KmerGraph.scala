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
}

/**
 * Graph showing connections between kmers.
 *
 * @param kmers k-mers which were observed.
 */
case class KmerGraph(protected val kmers: Array[Kmer]) {

  // source/sink kmers
  private def sourceKmers = kmers.filter(k => k.predecessors.length == 0)
  private def sinkKmers = kmers.filter(k => k.successors.length == 0)

  /**
   * Removes low frequency k-mers from the graph.
   */
  def removeLowFrequencyKmers(ratio: Double): KmerGraph = {
    val multiplicities = kmers.map(_.multiplicity).sorted()
    val threshold = (ratio *  multiplicities(multiplicities.length / 2).toDouble).toInt
    
    // the kmers to be filtered out
    val filterKmers = kmers.filter(_.multiplicity <= threshold)

    // map over kmers, and remove links to the kmers we are filtering
    filterKmers.foreach(kmer => {
      kmer.predecessors.foreach(_.removeLinks(kmer))
      kmer.successors.foreach(_.removeLinks(kmer))
    })

    // build new graph
    KmerGraph(kmers.filter(_.multiplicity > threshold))
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
