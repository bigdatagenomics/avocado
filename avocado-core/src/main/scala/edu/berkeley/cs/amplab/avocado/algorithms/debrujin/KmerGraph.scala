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

package edu.berkeley.cs.amplab.avocado.algorithms.debrujin

import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import scala.collection.mutable.{ArrayBuffer, HashMap, HashSet, PriorityQueue}

/**
 * For our purposes, a read is a list of kmers.
 *
 * @param record Read to build assembly read from.
 */
case class AssemblyRead (record: ADAMRecord) {
}

/**
 * A kmer prefix is a string of length k-1.
 *
 * @param string Kmer base string.
 */
case class KmerPrefix (string: String) {

  def equals(kp: KmerPrefix): Boolean = string == kp.string
}

/**
 * A kmer has a prefix of length k - 1 and a unit length suffix.
 *
 * @param prefix Prefix for kmer of length (k - 2)
 * @param suffix Last base of kmer
 */
case class Kmer (prefix: KmerPrefix, suffix: Char) {
  var reads = new HashSet[AssemblyRead]
  var mult: Int = 0
  var isCanon: Boolean = false

  def incrementMult() {
    mult += 1
  }

  def nextPrefix(): KmerPrefix = {
    KmerPrefix(prefix.string.drop(1) + suffix)
  }

  def equals(e: Kmer): Boolean = {
    prefix == e.prefix && suffix == e.suffix
  }

  override def toString(): String = {
    prefix.string + "[" + suffix + "]"
  }
}

/**
 *class representing a path made of kmers.
 *
 * @param edges Edges of kmer graph.
 */
class KmerPath (val edges: Seq[Kmer]) {

  val multSum: Int = edges.map(_.mult).fold(0)(_ + _)

  /**
   * Builds haplotype string from a kmer path using overlap between kmers.
   *
   * @return String representing haplotype from kmers.
   */
  def asHaplotypeString (): String = {
    var sb = ""
    if (edges.length > 0) {
      sb += edges(0).prefix.string
      for (e <- edges) {
        sb += e.suffix
      }
      sb
    } else {
      "" // TODO should throw exception
    }
  }

  def equals(kp: KmerPath): Boolean = {
    this.asHaplotypeString == kp.asHaplotypeString
  }
}

/**
 * An ordering representing how Kmer paths are ordered. Kmer paths are ordered
 * by increasing mult sum.
 */
object KmerPathOrdering extends Ordering[KmerPath] {

  /**
   * Compares two kmer paths, on the basis of mult sum.
   *
   * @param path1 Kmer path to evaluate.
   * @param path2 Kmer path to evaluate.
   * @return (-1, 0, 1) if path1 has mult sum (less than, equal, greater than) path2.
   */
  def compare (path1: KmerPath, path2: KmerPath): Int = {
    if (path1.multSum < path2.multSum) {
      -1
    } else if (path1.multSum > path2.multSum) {
      1
    } else {
      0
    }
  }
}

/**
 * Graph showing connections between kmers.
 *
 * @param kLen Kmer length.
 * @param readLen Read length.
 * @param regionLen Length of the active region.
 */
class KmerGraph (kLen: Int, readLen: Int, regionLen: Int, reference: String) {
  val spurThreshold = kLen // TODO(peter, 11/26) how to choose thresh?

  //var reads: HashMap[String,AssemblyRead] = null
  private var reads: Seq[AssemblyRead] = null

  // source/sink kmers
  val sourceKmer = new Kmer(new KmerPrefix(reference.take(kLen - 1)), reference.drop(kLen - 1).head)
  val sinkKmer = new Kmer(new KmerPrefix(reference.dropRight(1).takeRight(kLen - 1)),
                          reference.last)

  // The actual kmer graph consists of unique K-1 prefixes and kmers connected
  // by vertices.
  // TODO: Remove KmerPrefix class 
  var prefixes = new HashMap[String,KmerPrefix]
  var kmers = new HashMap[KmerPrefix,HashSet[Kmer]]

  // Paths through the kmer graph are in order of decreasing total mult.
  var allPaths = new PriorityQueue[KmerPath]()(KmerPathOrdering)

  /**
   * From a read, adds kmers to the graph.
   *
   * @param r Read to add kmers from.
   */
  def insertReadKmers (r: AssemblyRead): Unit = {
    // Construct L - K + 1 kmers, initially connected to the source and sink.
    val readSeq = r.record.getSequence.toString
    val offsets = 0 until readLen - kLen + 1
    val ks = offsets.map(idx => {
      val prefixStr = readSeq.substring(idx, idx + kLen - 1)
      val prefix = prefixes.getOrElseUpdate(prefixStr, new KmerPrefix(prefixStr))
      val suffix = readSeq.charAt(idx + kLen - 1)
      var k = new Kmer(prefix, suffix)

      // add kmer if seen
      val kmerSet = kmers.getOrElseUpdate(k.prefix, new HashSet[Kmer])
      if (kmerSet.forall(!_.equals(k))) {
        kmerSet += k
      }
      
      // increase multiplicity per kmer seen and attach read evidence
      kmerSet.filter(_.equals(k)).foreach(km => {
        km.incrementMult
        km.reads.add(r)
      })

      k
    })

    /*
    // Add a vertex in between each adjacent pair of kmers.
    (ks, ks.drop(1)).zipped.map((k1, k2) => {
      var v = new KmerVertex
      // TODO(peter, 11/26) check the order of k1 and k2!
      k1.right = v
      k2.left = v
      v.left.add(k1)
      v.right.add(k2)
    })

    // Add vertices at the ends.
    if (ks.head.equals(sourceKmer)) {
      ks.head.left = source
    } else {
      ks.head.left = new KmerVertex
    }
    ks.head.left.right.add(ks.head)
    if (ks.last.equals(sinkKmer)) {
      ks.last.right = sink
    } else {
      ks.last.right = new KmerVertex
    }
    ks.last.left.left.add(ks.last)
    */
  }

  /**
   * Inserts a group of reads into the kmer graph.
   *
   * @param readGroup Sequence of reads to add.
   */
  def insertReads (readGroup: Seq[ADAMRecord]): Unit = {
    reads = readGroup.map(x => new AssemblyRead(x))
    reads.foreach(r => insertReadKmers(r))
  }

  def exciseKmer (k: Kmer): Unit = {
    // TODO(peter, 11/27) for spur removal.
  }

  /**
   * Builds graph.
   */
  def connectGraph (): Unit = {
    /*
    // Consolidate equivalent kmers.
    for ((prefix, ks) <- kmers) {
      // Each equivalent (prefix, suffix) pair has an arbitrary "canonical" kmer.
      var canonKs = new HashMap[Char,Kmer]
      for (k <- ks) {
        var canonK = canonKs.getOrElseUpdate(k.suffix, k)
        canonK.isCanon = true
        if (k != canonK) {
          // Consolidate equivalent kmers together (i.e., same prefix and suffix)
          // with additive mults. Also fixup the vertices.
          mergeVertices(canonK.left, k.left)
          mergeVertices(canonK.right, k.right)
          exciseVertexKmer(canonK.left, k)
          exciseVertexKmer(canonK.right, k)
          canonK.reads ++= k.reads
          canonK.mult += k.mult
        }
      }
    }

    // Remove non-canonical kmers.
    for ((_, ks) <- kmers) {
      // TODO(peter, 12/5) if we keep references to non-canon kmers (e.g., in
      // AssemblyRead) we will have to get rid of those refs as well.
      ks --= ks.filter(k => !k.isCanon)
    }

    // Connect kmers to the source/sink when valid.
    for ((_, ks) <- kmers) {
      for (k <- ks) {
        if (k.left.left.size == 0) {
          k.left = source
        }
        if (k.right.right.size == 0) {
          k.right = sink
        }
      }
    }
    */
  }

  /**
   * Removes spurs from graph.
   */
  def removeSpurs (): Unit = {
    // Remove all kmer paths with length <= spurThreshold, connected to the
    // graph source and sink.
    // TODO(peter, 11/27) make sure the path lengths are initialized and
    // updated correctly!
    /*for (sk <- source.right) {
      var pathLen = 0
      var k = sk
      while (k.right.right.size == 1) {
        pathLen += 1
        k = k.right.right.head
      }
      if (pathLen <= spurThreshold) {
        k = k.left.left.head
        while (!source.equals(k)) {
          val lk = k.left.left.head
          exciseVertex(k.left)
          exciseKmer(k)
          k = lk
        }
      }
    }
    for (sk <- sink.left) {
      var pathLen = 0
      var k = sk
      while (k.left.left.size == 1) {
        pathLen += 1
        k = k.left.left.head
      }
      if (pathLen <= spurThreshold) {
        k = k.right.right.head
        while (!source.equals(k)) {
          val rk = k.right.right.head
          exciseVertex(k.right)
          exciseKmer(k)
          k = rk
        }
      }
    }*/
  }

  /**
   * Performs depth first search and scores all paths through the final graph.
   */
  def enumerateAllPaths (): Unit = {
    // TODO(peter, 12/9) arbitrary max assembly bound.
    //TODO change to use list cons
    val maxDepth = regionLen + readLen - kLen + 1
    var edges = new ArrayBuffer[Kmer]
    var paths = List[KmerPath]()

    def allPathsDFS (v: KmerPrefix, depth: Int): Unit = {
      if (v.equals(sinkKmer.nextPrefix)) {
        val path = new KmerPath(edges.clone.toSeq)
        paths = path :: paths
      } else if (depth <= maxDepth) {
        for (k <- kmers.getOrElse(v, Set[Kmer]())) {
          edges += k
          allPathsDFS(k.nextPrefix, depth + 1)
          edges.remove(edges.length - 1)
        }
      }
    }

    allPathsDFS(sourceKmer.prefix, 0)

    paths.foreach(p => allPaths.enqueue(p))
  }

  /**
   * Returns sorted priority queue containing all paths in the graph.
   *
   * @return Sorted priority queue of paths in the graph.
   */
  def getAllPaths (): PriorityQueue[KmerPath] = allPaths
}
