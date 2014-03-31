/*
 * Copyright (c) 2013. Regents of the University of California
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

package edu.berkeley.cs.amplab.avocado.calls.reads

import edu.berkeley.cs.amplab.adam.avro.{ADAMContig, 
                                         ADAMGenotype,
                                         ADAMGenotypeAllele,
                                         ADAMRecord, 
                                         ADAMVariant}
import edu.berkeley.cs.amplab.adam.models.ADAMVariantContext
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._
import edu.berkeley.cs.amplab.adam.rich.RichADAMRecord
import edu.berkeley.cs.amplab.adam.rich.RichADAMRecord._
import edu.berkeley.cs.amplab.adam.util.{MdTag, PhredUtils}
import edu.berkeley.cs.amplab.avocado.calls.VariantCallCompanion
import edu.berkeley.cs.amplab.avocado.stats.AvocadoConfigAndStats
import net.sf.samtools.{Cigar, CigarOperator, CigarElement, TextCigarCodec}
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.{RDD}
import scala.collection.JavaConversions._
import scala.collection.mutable.{ArrayBuffer, Buffer, HashMap, HashSet, PriorityQueue, StringBuilder}
import scala.math._

object VariantType extends scala.Enumeration {
  type VariantType = Value
  val SNP, MNP, Insertion, Deletion = Value
}

object ReadCallAssemblyPhaser extends VariantCallCompanion {

  val callName = "AssemblyPhaser"
  val debug = false

  def apply (stats: AvocadoConfigAndStats,
             config: SubnodeConfiguration): ReadCallAssemblyPhaser = {

    new ReadCallAssemblyPhaser()
  }
  
}

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

/**
 * Pairwise alignment HMM. See the Durbin textbook (1998), chapter 4.
 */
class HMMAligner {
  var refSequence: String = null
  var refAligned: String = null

  var testSequence: String = null
  var testAligned: String = null

  var testQualities: String = null

  var refOffset = 0

  var paddedRefLen = 0
  var paddedTestLen = 0

  var stride = 0
  var matSize = 1

  var eta = Double.NegativeInfinity

  // This uses the quick and dirty numbers from the Dindel (2011) paper.
  // TODO(peter, 12/7) I'm forgetting a factor of 3 somewhere...
  val mismatchPrior = -3.0 - log10(3.0)
  val matchPrior = log10(1.0 - 1.0e-3)
  val indelPrior = -4.0
  val indelToMatchPrior = -4.0
  val indelToIndelPrior = -4.0

  private var matches: Array[Double] = Array(0.0)
  private var inserts: Array[Double] = Array(0.0)
  private var deletes: Array[Double] = Array(0.0)

  private var traceMatches: Array[Char] = Array('\0')
  private var traceInserts: Array[Char] = Array('\0')
  private var traceDeletes: Array[Char] = Array('\0')

  private var alignment: ArrayBuffer[(Int, Char)] = null

  private var alignmentLikelihood = Double.NegativeInfinity
  private var alignmentPrior = Double.NegativeInfinity

  /**
   * Aligns sequences.
   *
   * @param refSequence Reference sequence over the active region.
   * @param testSequence Sequence being scored.
   * @param testQualities String of qualities. Not currently used.
   * @return True if sequence has variants.
   */
  def alignSequences (refSequence: String, testSequence: String, testQualities: String): Boolean = {

    def fpCompare(a: Double, b: Double, eps: Double): Boolean = abs(a - b) <= eps

    paddedRefLen = refSequence.length + 1
    paddedTestLen = testSequence.length + 1
    stride = paddedRefLen
    val oldMatSize = matSize
    matSize = max(matSize, paddedTestLen * paddedRefLen)
    if (matSize > oldMatSize && matSize > 0) {
      matches = new Array[Double](matSize)
      inserts = new Array[Double](matSize)
      deletes = new Array[Double](matSize)
      traceMatches = new Array[Char](matSize)
      traceInserts = new Array[Char](matSize)
      traceDeletes = new Array[Char](matSize)
    } else if (matSize <= 1 || oldMatSize <= 1) {
      matches = new Array[Double](1)
      inserts = new Array[Double](1)
      deletes = new Array[Double](1)
      traceMatches = new Array[Char](1)
      traceInserts = new Array[Char](1)
      traceDeletes = new Array[Char](1)
    }

    // Note: want to use the _test_ haplotype length here, not the ref length.
    eta = -log10(1.0 + testSequence.length)

    // Compute the optimal alignment.
    // TODO(peter, 12/4) shortcut b/w global and local alignment: use a custom
    // start position in the reference haplotype.
    matches(0) = 2.0 * eta
    inserts(0) = Double.NegativeInfinity
    deletes(0) = Double.NegativeInfinity
    for (i <- 0 until paddedTestLen) {
      for (j <- 0 until paddedRefLen) {
        if (i > 0 || j > 0) {
          val (m, trM) = if (i >= 1 && j >= 1) {
            val testBase = testSequence(i - 1)
            val refBase = refSequence(j - 1)
            // TODO(peter, 12/7) there is a second constant term to the prior...
            val prior = if (testBase == refBase) {
              matchPrior / (indelPrior * indelPrior)
            } else {
              mismatchPrior / (indelPrior * indelPrior)
            }
            val idx = (i - 1) * stride + (j - 1)
            val mMatch = matches(idx)
            val mInsert = inserts(idx)
            val mDelete = deletes(idx)
            val mMax = max(mMatch, max(mInsert, mDelete)) + prior
            val t = if (fpCompare(mMax, mMatch, 0.01)) {
              'M'
            } else if (fpCompare(mMax, mInsert, 0.01)) {
              'I'
            } else if (fpCompare(mMax, mDelete, 0.01)) {
              'D'
            } else {
              '.'
            }
            (mMax, t)
          } else {
            (Double.NegativeInfinity, '.')
          }
          val (ins, trIns) = if (i >= 1) {
            val idx = (i-1) * stride + j
            val insMatch = matches(idx) + indelToMatchPrior
            val insInsert = inserts(idx) + indelToIndelPrior
            val insMax = max(insMatch, insInsert)
            val t = if (fpCompare(insMax, insMatch, 0.01)) {
              'M'
            } else if (fpCompare(insMax, insInsert, 0.01)) {
              'I'
            } else {
              '.'
            }
            (insMax, t)
          } else {
            (Double.NegativeInfinity, '.')
          }
          val (del, trDel) = if (j >= 1) {
            val idx = i * stride + (j - 1)
            val delMatch = matches(idx) + indelToMatchPrior
            val delDelete = deletes(idx) + indelToIndelPrior
            val delMax = max(delMatch, delDelete)
            val t = if (fpCompare(delMax, delMatch, 0.01)) {
              'M'
            } else if (fpCompare(delMax, delDelete, 0.01)) {
              'D'
            } else {
              '.'
            }
            (delMax, t)
          } else {
            (Double.NegativeInfinity, '.')
          }
          val idx = i * stride + j
          matches(idx) = m
          inserts(idx) = ins
          deletes(idx) = del
          traceMatches(idx) = trM
          traceInserts(idx) = trIns
          traceDeletes(idx) = trDel
        }
      }
    }

    if (matSize > 0) {
      alignmentLikelihood = max(matches(matSize - 1), max(inserts(matSize - 1), deletes(matSize - 1)))
    } else {
      alignmentLikelihood = max(matches(0), max(inserts(0), deletes(0)))
    }

    def printArray(a: Array[Double]) {
      for (i <- 0 until paddedTestLen) {
        var s = ""
        for (j <- 0 until paddedRefLen) {
          val idx = i * stride + j
          s += "%2.2f" format a(idx)
          s += "\t"
        }
        if (ReadCallAssemblyPhaser.debug) println(s)
      }
    }
    
    if (ReadCallAssemblyPhaser.debug) {
      println("matches")
      printArray(matches)
      println("inserts")
      printArray(inserts)
      println("deletes")
      printArray(deletes)
    }

    // Traceback to get the aligned sequences.
    var hasVariants = false
    var numSnps = 0
    var numIndels = 0
    var revAlignment = ""
    var revAlignedTestSeq = ""
    var revAlignedRefSeq = ""

    def indexMax(a: (Double, Int), b: (Double, Int)): (Double, Int) = {
      if (a._1 > b._1) {
        a
      } else {
        b
      }
    }
    
    def getArrayMax(array: Array[Double]): (Double, Int) = {
      array.zipWithIndex.reduce(indexMax)
    }

    val highestM = getArrayMax(matches)
    val highestI = getArrayMax(inserts)
    val highestD = getArrayMax(deletes)
    val highestIdx = indexMax(highestM, indexMax(highestI, highestD))._2

    var i = paddedTestLen - 1
    var j = paddedRefLen - 1
    
    if (ReadCallAssemblyPhaser.debug)
      println(i + ", " + j)

    while (i > 0 && j > 0) {
      val idx = i * stride + j
      
      if (ReadCallAssemblyPhaser.debug) {
        println(i + ", " + j + (": %1.1f" format matches(idx)) + (", %1.1f" format inserts(idx)) +
                (", %1.1f" format deletes(idx)))
      }
              
      val bestScore = max(matches(idx), max(inserts(idx), deletes(idx)))
      // TODO(peter, 12/7) here, build the aligned sequences.
      val tr = if (fpCompare(bestScore, matches(idx), 0.0001)) {
        revAlignedTestSeq += testSequence(i - 1)
        revAlignedRefSeq += refSequence(j - 1)
        if (testSequence(i - 1) != refSequence(j - 1)) {
          hasVariants = true
          numSnps += 1
          revAlignment += 'X'
        }
        else {
          revAlignment += '='
        }
        // FIXME
        traceMatches(idx)
      } else if (fpCompare(bestScore, inserts(idx), 0.0001)) {
        revAlignedTestSeq += testSequence(i - 1)
        revAlignedRefSeq += '_'
        hasVariants = true
        numIndels += 1
        revAlignment += 'I'
        traceInserts(idx)
      } else if (fpCompare(bestScore, deletes(idx), 0.0001)) {
        revAlignedRefSeq += refSequence(j - 1)
        revAlignedTestSeq += '_'
        hasVariants = true
        numIndels += 1
        revAlignment += 'D'
        traceDeletes(idx)
      } else {
        revAlignedTestSeq += 'x'
        revAlignedRefSeq += 'x'
        '.'
      }

      tr match {
        case 'M' => {
          i -= 1
          j -= 1
        }
        case 'I' => {
          i -= 1
        }
        case 'D' => {
          j -= 1
        }
        case _ => {
          i -= 1
          j -= 1
        } // TODO(peter, 12/8) the alignment is bad (probably a bug).
      }
    }
    
    testAligned = revAlignedTestSeq.reverse
    refAligned = revAlignedRefSeq.reverse

    if (ReadCallAssemblyPhaser.debug) {
      println("ta: " + testAligned)
      println("ra: " + refAligned)
    }
      
    var unitAlignment = revAlignment.reverse
    if (ReadCallAssemblyPhaser.debug) println(unitAlignment)
    alignment = new ArrayBuffer[(Int, Char)]
    var alignSpan: Int = 0
    var alignMove: Char = '.'
    for (i <- 0 until unitAlignment.length) {
      val move = unitAlignment(i)
      if (move != alignMove) {
        if (alignSpan > 0) {
          val tok = (alignSpan, alignMove)
          alignment += tok
        }
        alignSpan = 1
        alignMove = move
      }
      else {
        alignSpan += 1
      }
    }
    val tok = (alignSpan, alignMove)
    alignment += tok

    // Compute the prior probability of the alignments, with the Dindel numbers.
    alignmentPrior = mismatchPrior * numSnps + indelPrior * numIndels
    if (ReadCallAssemblyPhaser.debug) println("ap: " + alignmentPrior)

    // Returns whether the alignment has any SNPs or indels.
    hasVariants
  }

  /**
   * Compute the (log10) likelihood of aligning the test sequence to the ref.
   *
   * @return Log10 likelihood of alignment.
   */
  def getLikelihood (): Double = alignmentLikelihood

  /**
   * Compute the (log10) prior prob of observing the given alignment.
   *
   * @return Prior for alignment.
   */
  def getPrior (): Double = alignmentPrior

  /**
   * Compute the alignment tokens (equivalent to cigar).
   *
   * @return Alignment tokens.
   */
  def getAlignment (): ArrayBuffer[(Int, Char)] = {
    alignment.clone
  }

  /**
   * Return the aligned sequences.
   *
   * @return Tuple of sequences aligned to (ref, test) sequences.
   */
  def getAlignedSequences (): (String, String) = (refAligned, testAligned)
}

/**
 * Haplotype generated from HMM alignment.
 *
 * @param sequence String representing haplotype alignment.
 */
class Haplotype (val sequence: String) {
  var perReadLikelihoods = new ArrayBuffer[Double]
  var readsLikelihood = Double.NegativeInfinity
  var hasVariants = false
  var alignment = new ArrayBuffer[(Int, Char)]
  
  /**
   * Score likelihood of reads when assembled into haplotype.
   *
   * @param hmm HMM aligner to use.
   * @param reads Sequence of reads to use.
   * @return Likelihood that reads are properly aligned.
   */
  def scoreReadsLikelihood (hmm: HMMAligner, reads: Seq[ADAMRecord]): Double = {
    perReadLikelihoods.clear
    readsLikelihood = 0.0
    for (r <- reads) {
      try {
        if (ReadCallAssemblyPhaser.debug) println(r.getSequence.toString + ", " + sequence)
        hmm.alignSequences(sequence, r.getSequence.toString, null)
        val readLike = hmm.getLikelihood // - hmm.getPriora
        perReadLikelihoods += readLike
        readsLikelihood += readLike
      } catch {
        case _ : Throwable => {
          perReadLikelihoods += 0.0
          readsLikelihood += 0.0
        }
      }
    }
    readsLikelihood
  }

  /**
   * Aligns reads to reference, and stores cigar details.
   *
   * @param hmm HMM aligner to use.
   * @param refHaplotype Haplotype for reference in this location.
   * @return True if region has variants.
   */
  def alignToReference (hmm: HMMAligner, refHaplotype: Haplotype): Boolean = {
    // TODO(peter, 12/8) store the alignment details (including the cigar).
    hasVariants = hmm.alignSequences(refHaplotype.sequence, sequence, null)
    alignment = hmm.getAlignment
    hasVariants
  }

  override def toString(): String = {
    sequence
  }
}

/**
 * Haplotypes are ordered by increasing reads likelihood, assuming they
 * come from the same group of reads.
 */
object HaplotypeOrdering extends Ordering[Haplotype] {

  /**
   * Compares two haplotypes. Returns (-1, 0, 1) if h1 has (lower, same, higher) read
   * likelihood than h2.
   *
   * @param h1 First haplotype to compare.
   * @param h2 Second haplotype to compare.
   * @return Ordering info for haplotypes.
   */
  def compare (h1: Haplotype, h2: Haplotype): Int = {
    if (h1.readsLikelihood < h2.readsLikelihood) {
      -1
    } else if (h1.readsLikelihood > h2.readsLikelihood) {
      1
    } else {
      0
    }
  }
}

/**
 * Class for a pairing of two haplotypes.
 *
 * @param haplotype1 First haplotype of pair.
 * @param haplotype2 Second haplotype of pair.
 */
class HaplotypePair (val haplotype1: Haplotype, val haplotype2: Haplotype) {
  var pairLikelihood = Double.NegativeInfinity
  var hasVariants = false

  override def toString(): String = {
    haplotype1.sequence + ", " + haplotype2.sequence + ", " + ("%1.3f" format pairLikelihood)
  }

  /**
   * Exponentiates two numbers by base of 10, adds together, takes base 10 log, and returns.
   *
   * @param x1 First digit to sum.
   * @param x2 Second digit to sum.
   * @return Sum of two digits after exponentation and logarithm.
   *
   * @see approxLogSumExp10
   */
  def exactLogSumExp10 (x1: Double, x2: Double): Double = {
    log10(pow(10.0, x1) + pow(10.0, x2))
  }

  /**
   * Exponentiates two numbers by base of 10, adds together, takes base 10 log, and returns.
   *
   * @param x1 First digit to sum.
   * @param x2 Second digit to sum.
   * @return Sum of two digits after exponentation and logarithm.
   *
   * @see exactLogSumExp10
   */
  def approxLogSumExp10 (x1: Double, x2: Double): Double = {
    exactLogSumExp10(x1, x2)
  }

  /**
   * Scores likelihood of two paired haplotypes and their alignment.
   *
   * @param hmm HMM aligner to use.
   * @param reads Sequence of reads that are evidence for haplotype.
   * @return Phred scaled likelihood.
   */
  def scorePairLikelihood (hmm: HMMAligner, reads: Seq[ADAMRecord]): Double = {
    var readsProb = 0.0
    for (i <- 0 until reads.length) {
      val readLikelihood1 = haplotype1.perReadLikelihoods(i)
      val readLikelihood2 = haplotype2.perReadLikelihoods(i)
      readsProb += exactLogSumExp10(readLikelihood1, readLikelihood2) - log10(2.0)
    }
    hmm.alignSequences(haplotype2.sequence, haplotype1.sequence, null)
    val priorProb = hmm.getPrior
    pairLikelihood = readsProb + priorProb
    pairLikelihood
  }

  /**
   * Aligns haplotype pair to reference. Joins variants of alignments.
   *
   * @param hmm HMM aligner to use.
   * @param refHaplotype Reference haplotype for active region.
   * @return True if region has variants.
   */
  def alignToReference (hmm: HMMAligner, refHaplotype: Haplotype): Boolean = {
    hasVariants = haplotype1.alignToReference(hmm, refHaplotype)
    if (haplotype2 != haplotype1) {
      hasVariants |= haplotype2.alignToReference(hmm, refHaplotype)
    }
    hasVariants
  }
}

/**
 * Haplotype pairs are ordered by increasing pairwise likelihood, assuming
 *  they come from the same read group.
 */
object HaplotypePairOrdering extends Ordering[HaplotypePair] {

  /**
   * Compares two haplotype pairs. Returns (-1, 0, 1) if first pair has (lower, same, higher)
   * pairwise likelihood.
   *
   * @param pair1 First haplotype pair to compare.
   * @param pair2 Second haplotype pair to compare.
   * @return Comparison of haplotype pairs.
   */
  def compare (pair1: HaplotypePair, pair2: HaplotypePair): Int = {
    if (pair1.pairLikelihood < pair2.pairLikelihood) {
      -1
    } else if (pair1.pairLikelihood > pair2.pairLikelihood) {
      1
    } else {
      0
    }
  }
}

/**
 * Phase (diploid) haplotypes with kmer assembly on active regions.
 */
class ReadCallAssemblyPhaser(val kmerLen: Int = 20,
                             val regionWindow: Int = 200) extends ReadCall {

  val companion = ReadCallAssemblyPhaser

  /**
   * From a read, returns the reference sequence.
   *
   * @param read Read from which to return sequence.
   * @return String containing reference sequence over this read.
   * 
   * @see https://github.com/bigdatagenomics/adam/blob/indel-realign/adam-commands/src/main/scala/edu/berkeley/cs/amplab/adam/util/MdTag.scala
   * @see getReference
   */
  def getReadReference (read: ADAMRecord): String = {
    val mdtag = MdTag(read.getMismatchingPositions.toString, read.getStart)

    val readSeq = RichADAMRecord(read).getSequence.toString
    val cigar = RichADAMRecord(read).samtoolsCigar

    var readPos = 0
    var refPos = 0
    var reference = ""

    val cigarEls: Buffer[CigarElement] = cigar.getCigarElements
    for (el <- cigarEls) {
      el.getOperator match {
        case CigarOperator.M => {
          for (i <- (0 until el.getLength)) {
            mdtag.mismatchedBase(refPos) match {
              case Some(b) => reference += b
              case None => reference += readSeq(readPos)
            }
            readPos += 1
            refPos += 1
          }
        }
        case CigarOperator.D => {
          for (i <- (0 until el.getLength)) {
            mdtag.deletedBase(refPos) match {
              case Some(b) => reference += b
              case None => {}
            }
            refPos += 1
          }
        }
        case CigarOperator.I => {
          readPos += el.getLength
        }
        case _ => {}
      }
    }

    reference
  }

  /**
   * Gets the reference from a set of reads. Works _provided_ that region has non-zero coverage
   * across whole region.
   *
   * @param region Sequence containing reads that cover region.
   * @return String of reference bases covering the region.
   *
   * @see getReadReference
   */
  def getReference (region: Seq[ADAMRecord]): String = {
    // TODO(peter, 12/5) currently, get the reference subsequence from the
    // MD tags of the ADAM records. Not entirely correct, because a region may
    // not be completely covered by reads, in which case the MD tags are
    // insufficient, so we ultimately want to resolve using the ref itself,
    // and not just the reads.
    val posRefs = (region.map(_.getStart), region.map(r => getReadReference(r)))
      .zipped.map((pos, ref) => (pos, ref))
      .sortBy(_._1)
    val startPos = posRefs(0)._1
    var reference = ""
    for ((pos, ref) <- posRefs) {
      // Here's an explanatory picture:
      //
      // OK:   [-----ref-----)
      //             [---read---)
      //
      // Skip: [-----ref-----)
      //         [---read---)
      //
      // Bail: [-----ref-----)
      //                         [---read---)
      val relPos = pos - startPos
      val offset = reference.length - relPos.toInt
      if (offset >= 0 && offset < ref.length) {
        try {
          reference += ref.substring(offset)
        } catch {
          case (e: StringIndexOutOfBoundsException) => {
            log.warn("String index out of bounds at: " + reference + ", " + ref + ", " + offset)
          }
        }
      } else if (offset < 0) {
        return ""
      }
    }
    reference
  }

  /**
   * Checks to see if region is active.
   *
   * @param region Sequence of reads over which to test.
   * @param ref Reference sequence over which to test.
   * @return True if region is active.
   */
  def isRegionActive (region: Seq[ADAMRecord], ref: String): Boolean = {
    // TODO(peter, 12/6) a very naive active region criterion. Upgrade asap!
    val activeLikelihoodThresh = -2.0
    var refHaplotype = new Haplotype(ref)
    var hmm = new HMMAligner
    val readsLikelihood = refHaplotype.scoreReadsLikelihood(hmm, region)
    readsLikelihood < activeLikelihoodThresh
  }

  /**
   * Performs assembly over a region.
   *
   * @param region Sequence of reads spanning the region.
   * @param ref String representing reference over the region.
   * @return Kmer graph corresponding to region.
   */
  def assemble (region: Seq[ADAMRecord], ref: String): KmerGraph = {
    val readLen = region(0).getSequence.length
    val regionLen = min(regionWindow + readLen - 1, ref.length)
    var kmerGraph = new KmerGraph(kmerLen, readLen, regionLen, ref)
    kmerGraph.insertReads(region)
    kmerGraph.connectGraph
    //kmer_graph.removeSpurs // TODO(peter, 11/27) debug: not doing spur removal atm.
    kmerGraph.enumerateAllPaths
    kmerGraph
  }

  /**
   * Emits a variant call, if a sequence is either heterozygous or homozygous non-ref.
   *
   * @param varType Type of variant.
   * @param varLength Length of variant.
   * @param varOffset Offset of variant in assembled haplotype.
   * @param refOffset Offset of variant versus reference haplotype.
   * @param varSequence Assembled sequence for variant.
   * @param refSequence Reference sequence.
   * @param heterozygousRef True if variant is heterozygous with a reference allele.
   * @param heterozygousNonref True if variant is heterozygous with two non-reference alleles.
   * @param phred Phred scaled variant quality.
   * @param sampleName Name of sample.
   * @param refName Name of reference sequence.
   * @param refId ID for reference.
   * @return List of genotypes.
   */
  def emitVariantCall (varType: VariantType.VariantType,
                       varLength: Int, 
                       varOffset: Int, 
                       refOffset: Int, 
                       varSequence: String, 
                       refSequence: String, 
                       heterozygousRef: Boolean, 
                       heterozygousNonref: Boolean, 
                       phred: Int, 
                       refPos: Long, 
                       sampleName: String, 
                       refName: String, 
                       refId: Int): List[ADAMGenotype] = {
    assert(!(heterozygousRef && heterozygousNonref))

    val refAllele = if (varType != VariantType.Insertion) {
      refSequence.substring(refOffset, refOffset + varLength)
    } else {
      ""
    }
    val altAllele = if (varType != VariantType.Deletion) {
      varSequence.substring(varOffset, varOffset + varLength)
    } else {
      ""
    }

    if (heterozygousRef) {
      val alleles = List(ADAMGenotypeAllele.Ref, ADAMGenotypeAllele.Alt)

      val contig = ADAMContig.newBuilder
        .setContigId(refId)
        .setContigName(refName)
        .build
      val variant = ADAMVariant.newBuilder
        .setContig(contig)
        .setReferenceAllele(refAllele)
        .setVariantAllele(altAllele)
        .setPosition(refPos + refOffset)
        .build
      val genotype = ADAMGenotype.newBuilder ()
        .setVariant(variant)
        .setSampleId (sampleName)
        .setGenotypeQuality(phred)
        .setExpectedAlleleDosage(1.0f)
        .setAlleles(alleles)
	.build()

      List(genotype)
    } else if (!heterozygousRef && !heterozygousNonref) {
      val alleles = List(ADAMGenotypeAllele.Alt, ADAMGenotypeAllele.Alt)

      val contig = ADAMContig.newBuilder
        .setContigId(refId)
        .setContigName(refName)
        .build
      val variant = ADAMVariant.newBuilder
        .setContig(contig)
        .setReferenceAllele(refAllele)
        .setVariantAllele(altAllele)
        .setPosition(refPos + refOffset)
        .build
      val genotype = ADAMGenotype.newBuilder ()
        .setVariant(variant)
        .setSampleId (sampleName)
        .setGenotypeQuality(phred)
        .setExpectedAlleleDosage(1.0f)
        .setAlleles(alleles)
	.build()

      List(genotype)
    } else {
      print("not calling")
      List[ADAMGenotype]()
    }
  }

  /**
   * Phasing the assembled haplotypes to call variants. See:
   *
   * C.A. Albers, G. Lunter, D.G. MacArthur, G. McVean, W.H. Ouwehand, R. Durbin.
   * "Dindel: Accurate indel calls from short-read data." Genome Research 21 (2011).
   *
   * @param region Sequence of reads covering region.
   * @param kmerGraph Graph of kmers to use for calling.
   * @param ref String for reference in the region.
   * @return List of variant contexts found in the region.
   */
  def phaseAssembly (region: Seq[ADAMRecord], 
                     kmerGraph: KmerGraph, 
                     ref: String): List[ADAMVariantContext] = {
    var refHaplotype = new Haplotype(ref)

    // Score all haplotypes against the reads.
    var hmm = new HMMAligner
    refHaplotype.scoreReadsLikelihood(hmm, region)
    var orderedHaplotypes = new PriorityQueue[Haplotype]()(HaplotypeOrdering)
    for (path <- kmerGraph.getAllPaths) {
      val haplotype = new Haplotype(path.asHaplotypeString)
      haplotype.scoreReadsLikelihood(hmm, region)
      orderedHaplotypes.enqueue(haplotype)
    }

    // Pick the top X-1 haplotypes and the reference haplotype.
    val maxNumBestHaplotypes = 16
    val numBestHaplotypes = min(maxNumBestHaplotypes, orderedHaplotypes.length)
    var bestHaplotypes = new ArrayBuffer[Haplotype]
    bestHaplotypes += refHaplotype
    for (i <- 1 to numBestHaplotypes) {
      bestHaplotypes += orderedHaplotypes.dequeue
    }

    // Score the haplotypes pairwise inclusively.
    var refHaplotypePair: HaplotypePair = null
    var orderedHaplotypePairs = new PriorityQueue[HaplotypePair]()(HaplotypePairOrdering)
    for (i <- 0 until bestHaplotypes.length) {
      for (j <- i until bestHaplotypes.length) {
        var pair = new HaplotypePair(bestHaplotypes(i), bestHaplotypes(j))
        pair.scorePairLikelihood(hmm, region)
        orderedHaplotypePairs.enqueue(pair)
        if (i == 0 && j == 0) {
          refHaplotypePair = pair
        }
      }
    }

    if (ReadCallAssemblyPhaser.debug) {
      println("After scoring, have:")
      orderedHaplotypePairs.foreach(println)
    }

    // Pick the best haplotype pairs with and without indels.
    val (calledHaplotypePair, uncalledHaplotypePair) = {
      var calledRes: HaplotypePair = null
      var uncalledRes: HaplotypePair = null
      do {
        val res = orderedHaplotypePairs.dequeue
        res.alignToReference(hmm, refHaplotype) match {
          case true => {
            if (calledRes == null) {
              calledRes = res
            }
          }
          case false => {
            if (uncalledRes == null) {
              uncalledRes = res
            }
          }
        }
      } while ((calledRes == null || uncalledRes == null) && orderedHaplotypePairs.length > 0)
      // TODO(peter, 12/8) this ought to be a pathological bug if it ever
      // happens (i.e., the ref-ref pair counts as having variants).
      // On the other hand, there might not be any valid variant haplotypes.
      // (FIXME: Really I should be using Option[].)
      if (uncalledRes == null) {
        uncalledRes = refHaplotypePair
      }
      (calledRes, uncalledRes)
    }

    // Compute the variant error probability and the equivalent phred score,
    // and use them for all called variants.
    val variantErrorProb = if (calledHaplotypePair != null) {
      val calledProbability = pow(10.0, calledHaplotypePair.pairLikelihood)
      val uncalledProbability = pow(10.0, uncalledHaplotypePair.pairLikelihood)
      uncalledProbability / (calledProbability + uncalledProbability)
    } else {
      1.0
    }
    val variantPhred = PhredUtils.successProbabilityToPhred(variantErrorProb)

    // TODO(peter, 12/8) Call variants, _without_ incorporating phasing for now.
    if (calledHaplotypePair != null) {
      var variants = List[ADAMGenotype]()
      val sortedRegion = region.sortBy(_.getStart)
      val firstRead = sortedRegion(0)
      val refPos = firstRead.getStart
      val refName = firstRead.getReferenceName.toString
      val refId = firstRead.getReferenceId
      val sampleName = firstRead.getRecordGroupSample.toString
      var calledHaplotypes = new HashSet[Haplotype]
      val calledHaplotype1 = calledHaplotypePair.haplotype1
      val calledHaplotype2 = calledHaplotypePair.haplotype2
      if (ReadCallAssemblyPhaser.debug) {
        println("Called:")
        println(calledHaplotype1 + ", " + calledHaplotype1.hasVariants + ", " + calledHaplotype1.alignment)
        println(calledHaplotype2 + ", " + calledHaplotype2.hasVariants + ", " + calledHaplotype2.alignment)
      }
      calledHaplotypes += calledHaplotypePair.haplotype1
      calledHaplotypes += calledHaplotypePair.haplotype2
      // FIXME this isn't quite right...
      val heterozygousRef = !calledHaplotypes.forall(_.hasVariants) && !calledHaplotypes.forall(!_.hasVariants)
      val heterozygousNonref = calledHaplotypes.filter(_.hasVariants).size > 1
      for (haplotype <- calledHaplotypes) {
        if (haplotype.hasVariants) {
          var variantOffset = 0
          var refOffset = 0
          for (tok <- haplotype.alignment) {
            val variantLength = tok._1
            val move = tok._2
            if (variantLength > 0 && (move == 'X' || move == 'I' || move == 'D')) {
              val variantType = move match {
                case 'X' => {
                  if (variantLength > 1) {
                    VariantType.MNP
                  } else {
                    VariantType.SNP
                  }
                }
                case 'I' => VariantType.Insertion
                case 'D' => VariantType.Deletion
              }
              val variant = emitVariantCall(variantType, variantLength,
                                            variantOffset, refOffset,
                                            haplotype.sequence, refHaplotype.sequence,
                                            heterozygousRef,
                                            heterozygousNonref,
                                            variantPhred,
                                            refPos,
                                            sampleName, refName, refId)
              variants = variants ::: variant
            }
            if (move != 'D') {
              variantOffset += variantLength
            }
            if (move != 'I') {
              refOffset += variantLength
            }
          }
        }
      }
      genotypesToVariantContext(variants)
    }
    else {
      List()
    }
  }

  /**
   * Method signature for variant calling operation.
   *
   * @param[in] pileupGroups An RDD containing reads.
   * @return An RDD containing called variants.
   */
  override def call (reads: RDD[ADAMRecord]): RDD[ADAMVariantContext] = {
    log.info("Grouping reads into active regions.")
    val activeRegions = reads.groupBy(r => r.getStart / regionWindow)
      .map(x => (getReference(x._2), x._2))
      .filter(x => isRegionActive(x._2, x._1))

    log.info("Calling variants with local assembly.")

    activeRegions.flatMap(x => {
      val ref = x._1
      val region = x._2
      if (ref.length > 0 && region.length > 0) {
        val kmerGraph = assemble(region, ref)
        phaseAssembly(region, kmerGraph, ref)
      }
      else {
        List()
      }
    })
  }

  override def isCallable (): Boolean = true
}
