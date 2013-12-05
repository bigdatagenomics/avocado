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

import edu.berkeley.cs.amplab.adam.avro.{ADAMRecord, ADAMVariant, ADAMGenotype, VariantType}
import edu.berkeley.cs.amplab.adam.util.{MdTag}
import edu.berkeley.cs.amplab.avocado.utils.{Phred}
import net.sf.samtools.{Cigar, CigarOperator, CigarElement, TextCigarCodec}
import org.apache.spark.{SparkContext, Logging}
import org.apache.spark.rdd.{RDD}
import scala.collection.JavaConversions._
import scala.collection.mutable.{ArrayBuffer, Buffer, HashMap, HashSet, PriorityQueue, StringBuilder}
import scala.math._

// For our purposes, a read is a list of kmers.
class AssemblyRead (_record: ADAMRecord) {
  val record = _record
  //var kmers: ArrayBuffer[Kmer] = null
}

// A kmer prefix is a string of length k-1.
class KmerPrefix(_string: String) {
  val string = _string
}

// A kmer has a prefix of length k-1 and a unit length suffix.
class Kmer (_prefix: KmerPrefix, _suffix: Char, _left: KmerVertex, _right: KmerVertex) {
  val prefix = _prefix
  val suffix = _suffix
  var left = _left
  var right = _right
  var reads = new HashSet[AssemblyRead]
  var mult: Int = 1
  var is_canon: Boolean = false
}

class KmerVertex {
  var left = new HashSet[Kmer]
  var right = new HashSet[Kmer]
}

class KmerPath (_edges: ArrayBuffer[Kmer]) {
  val edges = _edges
  var mult_sum: Int = 0
  for (e <- edges) {
    mult_sum += e.mult
  }

  def asHaplotypeString(): String = {
    var sb = new StringBuilder
    for (i <- Array.range(0, edges(0).prefix.string.length)) {
      sb += edges(0).prefix.string.charAt(i)
    }
    for (e <- edges) {
      sb += e.suffix
    }
    sb.toString
  }
}

object KmerPathOrdering extends Ordering[KmerPath] {
  // Kmer paths are ordered by increasing mult sum.
  def compare(path1: KmerPath, path2: KmerPath): Int = {
    if (path1.mult_sum < path2.mult_sum) {
      return -1
    }
    else if (path1.mult_sum > path2.mult_sum) {
      return 1
    }
    0
  }
}

// A mate pair is a pair of reads.
//class AssemblyMatePair (_r1: AssemblyRead, _r2: AssemblyRead) {
//  var read1: AssemblyRead = _r1
//  var read2: AssemblyRead = _r2
//}

class KmerGraph (_klen: Int, _read_len: Int, _region_len: Int) {
  val klen = _klen
  val read_len = _read_len
  val region_len = _region_len
  val spur_threshold = klen // TODO(peter, 11/26) how to choose thresh?

  //var reads: HashMap[String,AssemblyRead] = null
  private var reads: ArrayBuffer[AssemblyRead] = null

  // Convenient to explicitly have the graph source and sink.
  private val source = new KmerVertex
  private val sink = new KmerVertex

  // The actual kmer graph consists of unique K-1 prefixes and kmers connected
  // by vertices.
  private var prefixes = new HashMap[String,KmerPrefix]
  private var kmers = new HashMap[KmerPrefix,HashSet[Kmer]]

  // Paths through the kmer graph are in order of decreasing total mult.
  private var all_paths = new PriorityQueue[KmerPath]()(KmerPathOrdering)

  def insertReadKmers(r: AssemblyRead): Unit = {
    // Construct L-K+1 kmers, initially connected to the source and sink.
    val read_seq = r.record.getSequence.toString
    //val read_len = read_seq.length
    val offsets = ArrayBuffer.range(0, read_len - klen + 1)
    val ks = offsets.map(idx => {
      val prefix_str = read_seq.substring(idx, idx + klen - 1)
      val prefix = prefixes.getOrElseUpdate(prefix_str, new KmerPrefix(prefix_str))
      val suffix = read_seq.charAt(idx + klen - 1)
      var k = new Kmer(prefix, suffix, source, sink)
      k.reads.add(r)
      kmers.getOrElseUpdate(k.prefix, new HashSet[Kmer]) += k
      k
    })

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
    ks.head.left = new KmerVertex
    ks.head.left.right.add(ks.head)
    ks.last.right = new KmerVertex
    ks.last.left.left.add(ks.last)

    //r.kmers = ks
  }

  def insertReads(read_group: Seq[ADAMRecord]): Unit = {
    reads = ArrayBuffer(read_group.map(x => new AssemblyRead(x)) : _*)
    reads.map(r => insertReadKmers(r))
  }

  def mergeVertices(v1: KmerVertex, v2: KmerVertex): Unit = {
    // Merge v2 into v1.
    v2.left.map(k => v1.left.add(k))
    v2.right.map(k => v1.right.add(k))
  }

  def exciseVertexKmer(v: KmerVertex, k: Kmer): Unit = {
    v.left.remove(k)
    v.right.remove(k)
  }

  def exciseKmer(k: Kmer): Unit = {
    // TODO(peter, 11/27) for spur removal.
  }

  def exciseVertex(v: KmerVertex): Unit = {
    // TODO(peter, 11/27) for spur removal.
  }

  def connectGraph(): Unit = {
    // Consolidate equivalent kmers.
    for ((prefix, ks) <- kmers) {
      // Each equivalent (prefix, suffix) pair has an arbitrary "canonical" kmer.
      var canon_ks = new HashMap[Char,Kmer]
      for (k <- ks) {
        var canon_k = canon_ks.getOrElseUpdate(k.suffix, k)
        canon_k.is_canon = true
        if (k != canon_k) {
          // Consolidate equivalent kmers together (i.e., same prefix and suffix)
          // with additive mults. Also fixup the vertices.
          mergeVertices(canon_k.left, k.left)
          mergeVertices(canon_k.right, k.right)
          exciseVertexKmer(canon_k.left, k)
          exciseVertexKmer(canon_k.right, k)
          canon_k.reads ++= k.reads
          canon_k.mult += k.mult
        }
      }
    }

    // Remove non-canonical kmers.
    for ((_, ks) <- kmers) {
      // TODO(peter, 12/5) if we keep references to non-canon kmers (e.g., in
      // AssemblyRead) we will have to get rid of those refs as well.
      ks --= ks.filter(k => !k.is_canon)
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
  }

  def removeSpurs(): Unit = {
    // Remove all kmer paths with length <= spur_threshold, connected to the
    // graph source and sink.
    // TODO(peter, 11/27) make sure the path lengths are initialized and
    // updated correctly!
    for (sk <- source.right) {
      var path_len = 0
      var k = sk
      while (k.right.right.size == 1) {
        path_len += 1 //k.suffix.length
        k = k.right.right.head
      }
      if (path_len <= spur_threshold) {
        k = k.left.left.head
        while (k != source) {
          val lk = k.left.left.head
          exciseVertex(k.left)
          exciseKmer(k)
          k = lk
        }
      }
    }
    for (sk <- sink.left) {
      var path_len = 0
      var k = sk
      while (k.left.left.size == 1) {
        path_len += 1 //k.left.left.head.suffix.length
        k = k.left.left.head
      }
      if (path_len <= spur_threshold) {
        k = k.right.right.head
        while (k != source) {
          val rk = k.right.right.head
          exciseVertex(k.right)
          exciseKmer(k)
          k = rk
        }
      }
    }
  }

  //def threadMatePairs(): Unit = {
  //}

  def enumerateAllPaths(): Unit = {
    // Do DFS to enumerate and score all paths through the final graph.
    // TODO(peter, 12/9) arbitrary max assembly bound.
    val max_depth = region_len + read_len - klen + 1
    var edges = new ArrayBuffer[Kmer]
    def allPathsDFS(v: KmerVertex, depth: Int): Unit = {
      if (v == sink) {
        val path = new KmerPath(edges)
        all_paths.enqueue(path)
      }
      else if (depth <= max_depth) {
        for (k <- v.right) {
          edges += k
          allPathsDFS(k.right, depth + 1)
          edges.remove(edges.length - 1)
        }
      }
    }
    allPathsDFS(source, 0)
  }

  def getAllPaths(): PriorityQueue[KmerPath] = all_paths
}

/**
 * Pairwise alignment HMM. See the Durbin textbook (1998), chapter 4.
 */
class HMMAligner {
  var ref_sequence: String = null
  var ref_aligned: String = null

  var test_sequence: String = null
  var test_aligned: String = null

  var test_qualities: String = null

  var ref_offset = 0

  var padded_ref_len = 0
  var padded_test_len = 0

  var stride = 0
  var mat_size = 0

  var eta = Double.NegativeInfinity

  // This uses the quick and dirty numbers from the Dindel (2011) paper.
  // TODO(peter, 12/7) I'm forgetting a factor of 3 somewhere...
  val mismatch_prior = -3.0 - log10(3.0)
  val match_prior = log10(1.0 - 1.0e-3)
  val indel_prior = -4.0
  val indel_to_match_prior = -4.0
  val indel_to_indel_prior = -4.0

  private var matches: ArrayBuffer[Double] = null
  private var inserts: ArrayBuffer[Double] = null
  private var deletes: ArrayBuffer[Double] = null

  private var trace_matches: ArrayBuffer[Char] = null
  private var trace_inserts: ArrayBuffer[Char] = null
  private var trace_deletes: ArrayBuffer[Char] = null

  private var alignment: ArrayBuffer[(Int, Char)] = null

  private var alignment_likelihood = Double.NegativeInfinity
  private var alignment_prior = Double.NegativeInfinity

  //def alignSequences(_ref: String, _test: String, _ref_offset: Int): Unit = {
  def alignSequences(_ref: String, _test: String, _test_quals: String): Boolean = {
    ref_sequence = _ref
    test_sequence = _test
    test_qualities = _test_quals

    //ref_offset = _ref_offset

    padded_ref_len = ref_sequence.length + 1
    padded_test_len = test_sequence.length + 1
    stride = padded_ref_len
    val old_mat_size = mat_size
    mat_size = max(mat_size, padded_test_len * padded_ref_len)
    if (mat_size > old_mat_size) {
      matches = new ArrayBuffer[Double](mat_size)
      inserts = new ArrayBuffer[Double](mat_size)
      deletes = new ArrayBuffer[Double](mat_size)
      trace_matches = new ArrayBuffer[Char](mat_size)
      trace_inserts = new ArrayBuffer[Char](mat_size)
      trace_deletes = new ArrayBuffer[Char](mat_size)
    }

    // Note: want to use the _test_ haplotype length here, not the ref length.
    eta = -log10(1.0 + test_sequence.length)

    // Compute the optimal alignment.
    // TODO(peter, 12/4) shortcut b/w global and local alignment: use a custom
    // start position in the reference haplotype.
    matches(0) = 2.0 * eta
    inserts(0) = Double.NegativeInfinity
    deletes(0) = Double.NegativeInfinity
    for (i <- 0 to padded_test_len) {
      //for (j <- 0 to (padded_ref_len - ref_offset)) {
      for (j <- 0 to padded_ref_len) {
        if (i > 0 || j > 0) {
          val (m, tr_m) = if (i >= 1 && j >= 1) {
            val test_base = test_sequence(i-1)
            val ref_base = ref_sequence(j-1)
            // TODO(peter, 12/7) there is a second constant term to the prior...
            val prior = if (test_base == ref_base) {
              match_prior / (indel_prior * indel_prior)
            } else {
              mismatch_prior / (indel_prior * indel_prior)
            }
            val idx = (i-1) * stride + (j-1)
            val m_match = matches(idx)
            val m_insert = inserts(idx)
            val m_delete = deletes(idx)
            val m = max(m_match, max(m_insert, m_delete)) + prior
            val t = if (m == m_match) {
              'M'
            } else if (m == m_insert) {
              'I'
            } else if (m == m_delete) {
              'D'
            } else {
              '.'
            }
            (m, t)
          } else {
            (Double.NegativeInfinity, '.')
          }
          val (ins, tr_ins) = if (i >= 1) {
            val idx = (i-1) * stride + j
            val ins_match = matches(idx) + indel_to_match_prior
            val ins_insert = inserts(idx) + indel_to_indel_prior
            val ins = max(ins_match, ins_insert)
            val t = if (ins == ins_match) {
              'M'
            } else if (ins == ins_insert) {
              'I'
            } else {
              '.'
            }
            (ins, t)
          } else {
            (Double.NegativeInfinity, '.')
          }
          val (del, tr_del) = if (j >= 1) {
            val idx = i * stride + (j-1)
            val del_match = matches(idx) + indel_to_match_prior
            val del_delete = deletes(idx) + indel_to_indel_prior
            val del = max(del_match, del_delete)
            val t = if (del == del_match) {
              'M'
            } else if (del == del_delete) {
              'D'
            } else {
              '.'
            }
            (del, t)
          } else {
            (Double.NegativeInfinity, '.')
          }
          val idx = i * stride + j
          matches(idx) = m
          inserts(idx) = ins
          deletes(idx) = del
          trace_matches(idx) = tr_m
          trace_inserts(idx) = tr_ins
          trace_deletes(idx) = tr_del
        }
      }
    }
    alignment_likelihood = max(matches(mat_size - 1), max(inserts(mat_size - 1), deletes(mat_size - 1)))

    // Traceback to get the aligned sequences.
    var has_variants = false
    var num_snps = 0
    var num_indels = 0
    var rev_alignment = ""
    var rev_aligned_test_seq = ""
    var rev_aligned_ref_seq = ""
    var i = padded_test_len - 1
    var j = padded_ref_len - 1
    while (i > 0 && j > 0) {
      val idx = i * stride + j
      val best_score = max(matches(idx), max(inserts(idx), deletes(idx)))
      // TODO(peter, 12/7) here, build the aligned sequences.
      val tr = if (best_score == matches(idx)) {
        rev_aligned_test_seq += test_sequence(i-1)
        rev_aligned_ref_seq += ref_sequence(j-1)
        if (test_sequence(i-1) != ref_sequence(j-1)) {
          has_variants = true
          num_snps += 1
          rev_alignment += 'X'
        }
        else {
          rev_alignment += '='
        }
        trace_matches(idx)
      } else if (best_score == inserts(idx)) {
        rev_aligned_test_seq += test_sequence(i-1)
        has_variants = true
        num_indels += 1
        rev_alignment += 'I'
        trace_inserts(idx)
      } else if (best_score == deletes(idx)) {
        rev_aligned_ref_seq += ref_sequence(j-1)
        has_variants = true
        num_indels += 1
        rev_alignment += 'D'
        trace_deletes(idx)
      } else {
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
        case _ => {} // TODO(peter, 12/8) the alignment is bad (probably a bug).
      }
    }
    test_aligned = rev_aligned_test_seq.reverse
    ref_aligned = rev_aligned_ref_seq.reverse
    var unit_alignment = rev_alignment.reverse
    alignment = new ArrayBuffer[(Int, Char)]
    var align_span: Int = 0
    var align_move: Char = '.'
    for (i <- 0 to unit_alignment.length) {
      val move = unit_alignment(i)
      if (move != align_move) {
        if (align_span > 0) {
          val tok = (align_span, align_move)
          alignment += tok
        }
        align_span = 1
        align_move = move
      }
      else {
        align_span += 1
      }
    }
    val tok = (align_span, align_move)
    alignment += tok

    // Compute the prior probability of the alignments, with the Dindel numbers.
    alignment_prior = mismatch_prior * num_snps + indel_prior * num_indels

    // Returns whether the alignment has any SNPs or indels.
    has_variants
  }

  // Compute the (log10) likelihood of aligning the test sequence to the ref.
  def getLikelihood(): Double = alignment_likelihood

  // Compute the (log10) prior prob of observing the given alignment.
  def getPrior(): Double = alignment_prior

  // Compute the alignment tokens (equivalent to cigar).
  def getAlignment(): ArrayBuffer[(Int, Char)] = {
    alignment.clone
  }

  // Return the aligned sequences.
  def getAlignedSequences(): (String, String) = (ref_aligned, test_aligned)
}

class Haplotype (_sequence: String) {
  val sequence = _sequence
  var per_read_likelihoods = new ArrayBuffer[Double]
  var reads_likelihood = Double.NegativeInfinity
  var has_variants = false
  var alignment = new ArrayBuffer[(Int, Char)]

  def scoreReadsLikelihood(hmm: HMMAligner, reads: Seq[ADAMRecord]): Double = {
    per_read_likelihoods.clear
    reads_likelihood = 0.0
    for (r <- reads) {
      hmm.alignSequences(sequence, r.getSequence.toString, null)
      val read_like = hmm.getLikelihood // - hmm.getPriora
      per_read_likelihoods += read_like
      reads_likelihood += read_like
    }
    reads_likelihood
  }

  def alignToReference(hmm: HMMAligner, ref_haplotype: Haplotype): Boolean = {
    // TODO(peter, 12/8) store the alignment details (including the cigar).
    has_variants = hmm.alignSequences(ref_haplotype.sequence, sequence, null)
    alignment = hmm.getAlignment
    has_variants
  }
}

object HaplotypeOrdering extends Ordering[Haplotype] {
  // Haplotypes are ordered by increasing reads likelihood, assuming they
  // come from the same read group.
  def compare(h1: Haplotype, h2: Haplotype): Int = {
    if (h1.reads_likelihood < h2.reads_likelihood) {
      return -1
    }
    else if (h1.reads_likelihood > h2.reads_likelihood) {
      return 1
    }
    0
  }
}

class HaplotypePair (_h1: Haplotype, _h2: Haplotype) {
  val haplotype1 = _h1
  val haplotype2 = _h2
  var pair_likelihood = Double.NegativeInfinity
  var has_variants = false

  def exactLogSumExp10(x1: Double, x2: Double): Double = {
    log10(pow(10.0, x1) + pow(10.0, x2))
  }

  def approxLogSumExp10(x1: Double, x2: Double): Double = {
    // TODO(peter, 12/7)
    log10(pow(10.0, x1) + pow(10.0, x2))
  }

  def scorePairLikelihood(hmm: HMMAligner, reads: Seq[ADAMRecord]): Double = {
    var reads_prob = 0.0
    for (i <- 0 to reads.length) {
      val read_likelihood1 = haplotype1.per_read_likelihoods(i)
      val read_likelihood2 = haplotype2.per_read_likelihoods(i)
      reads_prob += exactLogSumExp10(read_likelihood1, read_likelihood2) - log10(2.0)
    }
    hmm.alignSequences(haplotype2.sequence, haplotype1.sequence, null)
    val prior_prob = hmm.getPrior
    pair_likelihood = reads_prob + prior_prob
    pair_likelihood
  }

  def alignToReference(hmm: HMMAligner, ref_haplotype: Haplotype): Boolean = {
    has_variants = haplotype1.alignToReference(hmm, ref_haplotype)
    if (haplotype2 != haplotype1) {
      has_variants |= haplotype2.alignToReference(hmm, ref_haplotype)
    }
    has_variants
  }
}

object HaplotypePairOrdering extends Ordering[HaplotypePair] {
  // Haplotype pairs are ordered by increasing pairwise likelihood, assuming
  // they come from the same read group.
  def compare(pair1: HaplotypePair, pair2: HaplotypePair): Int = {
    if (pair1.pair_likelihood < pair2.pair_likelihood) {
      return -1
    }
    else if (pair1.pair_likelihood > pair2.pair_likelihood) {
      return 1
    }
    0
  }
}

/**
 * Phase (diploid) haplotypes with kmer assembly on active regions.
 */
class ReadCallAssemblyPhaser extends ReadCall {
  override val callName = "AssemblyPhaser"

  val kmer_len = 20
  val region_window = 200

  val cigar_codec = new TextCigarCodec

  // See: (https://github.com/bigdatagenomics/adam/blob/indel-realign/adam-commands/src/main/scala/edu/berkeley/cs/amplab/adam/util/MdTag.scala).
  def getReadReference(read: ADAMRecord): String = {
    val mdtag = MdTag(read.getMismatchingPositions.toString, read.getStart)

    val read_seq = read.getSequence.toString
    val cigar = cigar_codec.decode(read.getCigar.toString)

    var read_pos = 0
    var ref_pos = 0
    var reference = ""

    val cigar_els: Buffer[CigarElement] = cigar.getCigarElements
    for (el <- cigar_els) {
      el.getOperator match {
        case CigarOperator.M => {
          for (i <- (0 until el.getLength)) {
            mdtag.mismatchedBase(ref_pos) match {
              case Some(b) => reference += b
              case None => reference += read_seq(read_pos)
            }
            read_pos += 1
            ref_pos += 1
          }
        }
        case CigarOperator.D => {
          for (i <- (0 until el.getLength)) {
            mdtag.deletedBase(ref_pos) match {
              case Some(b) => reference += b
              case None => {}
            }
            ref_pos += 1
          }
        }
        case CigarOperator.I => {
          read_pos += el.getLength
        }
        case _ => {}
      }
    }

    reference
  }

  def getReference(region: Seq[ADAMRecord]): String = {
    // TODO(peter, 12/5) currently, get the reference subsequence from the
    // MD tags of the ADAM records. Not entirely correct, because a region may
    // not be completely covered by reads, in which case the MD tags are
    // insufficient, so we ultimately want to resolve using the ref itself,
    // and not just the reads.
    val pos_refs = (region.map(_.getStart), region.map(r => getReadReference(r)))
                   .zipped.map((pos, ref) => (pos, ref))
                   .sortBy(_._1)
    val start_pos = pos_refs(0)._1
    var reference = ""
    for ((pos, ref) <- pos_refs) {
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
      val rel_pos = pos - start_pos
      val offset = reference.length - rel_pos.toInt
      if (offset >= 0) {
        reference += ref.substring(offset, ref.length)
      }
      else if (offset < 0) {
        return ""
      }
    }
    reference
  }

  def regionIsActive(region: Seq[ADAMRecord], ref: String): Boolean = {
    // TODO(peter, 12/6) a very naive active region criterion. Upgrade asap!
    val active_likelihood_thresh = -2.0
    var ref_haplotype = new Haplotype(ref)
    var hmm = new HMMAligner
    val reads_likelihood = ref_haplotype.scoreReadsLikelihood(hmm, region)
    reads_likelihood < active_likelihood_thresh
  }

  def assemble(region: Seq[ADAMRecord], ref: String): KmerGraph = {
    val read_len = region(0).getSequence.length
    val region_len = min(region_window + read_len - 1, ref.length)
    var kmer_graph = new KmerGraph(kmer_len, read_len, region_len)
    kmer_graph.insertReads(region)
    kmer_graph.connectGraph
    //kmer_graph.removeSpurs // TODO(peter, 11/27) debug: not doing spur removal atm.
    kmer_graph.enumerateAllPaths
    kmer_graph
  }

  def emitVariantCall(var_type: VariantType, var_length: Int, var_offset: Int, ref_offset: Int, var_sequence: String, ref_sequence: String, heterozygous_ref: Boolean, heterozygous_nonref: Boolean, phred: Int, ref_pos: Long, sample_name: String, ref_name: String): (ADAMVariant, List[ADAMGenotype]) = {
    val (genotype_str, allele_count) = if (heterozygous_ref) {
      ("0,1", 2)
    } else if (!heterozygous_ref && !heterozygous_nonref) {
      ("1,1", 2)
    } else {
      ("1,2", 3)
    }
    val ref_allele = if (var_type != VariantType.Insertion) {
      ref_sequence.substring(ref_offset, ref_offset + var_length)
    } else {
      ""
    }
    val alt_allele = if (var_type != VariantType.Deletion) {
      var_sequence.substring(var_offset, var_offset + var_length)
    } else {
      ""
    }
    val genotype = ADAMGenotype.newBuilder
        .setSampleId(sample_name)
        .setGenotype(genotype_str)
        .setPhredLikelihoods(phred.toString)
        .build
    val variant = ADAMVariant.newBuilder
        .setReferenceName(ref_name)
        .setStartPosition(ref_pos + ref_offset)
        .setReferenceAllele(ref_allele)
        .setAlternateAlleles(alt_allele)
        .setAlleleCount(allele_count)
        .setType(var_type)
        .build
    (variant, List(genotype))
  }

  /**
   * Phasing the assembled haplotypes to call variants. See:
   *
   * C.A. Albers, G. Lunter, D.G. MacArthur, G. McVean, W.H. Ouwehand, R. Durbin.
   * "Dindel: Accurate indel calls from short-read data." Genome Research 21 (2011).
   */
  def phaseAssembly(region: Seq[ADAMRecord], kmer_graph: KmerGraph, ref: String): Seq[(ADAMVariant, List[ADAMGenotype])] = {
    var ref_haplotype = new Haplotype(ref)

    // Score all haplotypes against the reads.
    var hmm = new HMMAligner
    ref_haplotype.scoreReadsLikelihood(hmm, region)
    var ordered_haplotypes = new PriorityQueue[Haplotype]()(HaplotypeOrdering)
    for (path <- kmer_graph.getAllPaths) {
      val haplotype = new Haplotype(path.asHaplotypeString)
      haplotype.scoreReadsLikelihood(hmm, region)
      ordered_haplotypes.enqueue(haplotype)
    }

    // Pick the top X-1 haplotypes and the reference haplotype.
    val max_num_best_haplotypes = 16
    val num_best_haplotypes = min(max_num_best_haplotypes, ordered_haplotypes.length)
    var best_haplotypes = new ArrayBuffer[Haplotype]
    best_haplotypes += ref_haplotype
    for (i <- 1 to num_best_haplotypes) {
      best_haplotypes += ordered_haplotypes.dequeue
    }

    // Score the haplotypes pairwise inclusively.
    var ref_haplotype_pair: HaplotypePair = null
    var ordered_haplotype_pairs = new PriorityQueue[HaplotypePair]()(HaplotypePairOrdering)
    for (i <- 0 to best_haplotypes.length) {
      for (j <- 0 to (i + 1)) {
        var pair = new HaplotypePair(best_haplotypes(i), best_haplotypes(j))
        pair.scorePairLikelihood(hmm, region)
        ordered_haplotype_pairs.enqueue(pair)
        if (i == 0 && j == 0) {
          ref_haplotype_pair = pair
        }
      }
    }

    // Pick the best haplotype pairs with and without indels.
    val (called_haplotype_pair, uncalled_haplotype_pair) = {
      var called_res: HaplotypePair = null
      var uncalled_res: HaplotypePair = null
      do {
        val res = ordered_haplotype_pairs.dequeue
        res.alignToReference(hmm, ref_haplotype) match {
          case true => {
            if (called_res == null) {
              called_res = res
            }
          }
          case false => {
            if (uncalled_res == null) {
              uncalled_res = res
            }
          }
        }
      } while ((called_res == null || uncalled_res == null) && ordered_haplotype_pairs.length > 0)
      // TODO(peter, 12/8) this ought to be a pathological bug if it ever
      // happens (i.e., the ref-ref pair counts as having variants).
      // On the other hand, there might not be any valid variant haplotypes.
      // (FIXME: Really I should be using Option[].)
      if (uncalled_res == null) {
        uncalled_res = ref_haplotype_pair
      }
      (called_res, uncalled_res)
    }

    // Compute the variant error probability and the equivalent phred score,
    // and use them for all called variants.
    val variant_error_prob = if (called_haplotype_pair != null) {
      val called_probability = pow(10.0, called_haplotype_pair.pair_likelihood)
      val uncalled_probability = pow(10.0, uncalled_haplotype_pair.pair_likelihood)
      uncalled_probability / (called_probability + uncalled_probability)
    } else {
      1.0
    }
    val variant_phred = Phred.probabilityToPhred(variant_error_prob)

    // TODO(peter, 12/8) Call variants, _without_ incorporating phasing for now.
    if (called_haplotype_pair != null) {
      var variants = new ArrayBuffer[(ADAMVariant, List[ADAMGenotype])]
      val sorted_region = region.sortBy(_.getStart)
      val first_read = sorted_region(0)
      val ref_pos = first_read.getStart
      val ref_name = first_read.getReferenceName.toString
      val sample_name = first_read.getRecordGroupSample.toString
      var called_haplotypes = new HashSet[Haplotype]
      val called_haplotype1 = called_haplotype_pair.haplotype1
      val called_haplotype2 = called_haplotype_pair.haplotype2
      called_haplotypes += called_haplotype_pair.haplotype1
      called_haplotypes += called_haplotype_pair.haplotype2
      val heterozygous_ref = called_haplotypes.map(!_.has_variants).foldLeft(false)((z, x) => z | x)
      val heterozygous_nonref = called_haplotypes.size > 1
      for (haplotype <- called_haplotypes) {
        if (haplotype.has_variants) {
          var variant_offset = 0
          var ref_offset = 0
          for (tok <- haplotype.alignment) {
            val variant_length = tok._1
            val move = tok._2
            if (variant_length > 0 && (move == 'X' || move == 'I' || move == 'D')) {
              val variant_type = move match {
                case 'X' => {
                  if (variant_length > 1) {
                    VariantType.MNP
                  } else {
                    VariantType.SNP
                  }
                }
                case 'I' => VariantType.Insertion
                case 'D' => VariantType.Deletion
              }
              val variant = emitVariantCall(variant_type, variant_length,
                                            variant_offset, ref_offset,
                                            haplotype.sequence, ref_haplotype.sequence,
                                            heterozygous_ref,
                                            heterozygous_nonref,
                                            variant_phred,
                                            ref_pos,
                                            sample_name, ref_name)
              variants += variant
            }
            if (move != 'D') {
              variant_offset += variant_length
            }
            if (move != 'I') {
              ref_offset += variant_length
            }
          }
        }
      }
      variants
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
  override def call(reads: RDD[ADAMRecord]): RDD[(ADAMVariant, List[ADAMGenotype])] = {
    log.info("Grouping reads into active regions.")
    val active_regions = reads.groupBy(r => r.getStart / region_window)
                              .map(x => (getReference(x._2), x._2))
                              .filter(x => regionIsActive(x._2, x._1))
    log.info("Found " + active_regions.count.toString + " active regions.")
    log.info("Calling variants with local assembly.")
    active_regions.flatMap(x => {
      val ref = x._1
      val region = x._2
      if (ref.length > 0 && region.length > 0) {
        val kmer_graph = assemble(region, ref)
        phaseAssembly(region, kmer_graph, ref)
      }
      else {
        List()
      }
    })
  }

  override def isCallable(): Boolean = true
}
