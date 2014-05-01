/*
 * Copyright (c) 2014. Mount Sinai School of Medicine
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use file except in compliance with the License.
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

/**
 * A kmer has a prefix of length k - 1 and a unit length suffix.
 *
 * @param kmerSeq sequence that will be split into k-1 len prefix and single char suffix
 * @param weight
 */
case class Kmer(val kmerSeq: String, val weight: Int = 1) {

  val k = kmerSeq.size
  val prefix: String = kmerSeq.substring(0, k - 1)
  val suffix: Char = kmerSeq.charAt(k - 1)

  def nextPrefix: String = {
    prefix.drop(1) + suffix
  }

  override def toString: String = {
    prefix + "[" + suffix + "]"
  }

  /**
   * Prints this kmer's connectivity in graphviz format. This is used by the KmerGraph class
   * for printing the graph as a graphviz graph.
   *
   * @return Returns a string describing the directed connectivity between this _k_-mer and
   * the next _k_-mer it points at.
   */
  def toDot: String = {
    prefix + " -> " + nextPrefix + " ;"
  }
}

/**
 * class representing a path made of kmers.
 *
 * @param edges Edges of kmer graph.
 */
class KmerPath(val edges: Seq[Kmer]) extends Ordered[KmerPath] {

  val weight: Int = edges.map(_.weight).reduce(_ + _)
  val len = edges.size
  /**
   * Builds haplotype string from a kmer path using overlap between kmers.
   *
   * @return String representing haplotype from kmers.
   */
  lazy val haplotypeString: String = {
    val builder = new StringBuilder
    if (len > 0) {
      builder.append(edges(0).prefix)
      edges.foreach(edge => builder.append(edge.suffix))
    }
    builder.toString
  }

  def equals(kp: KmerPath): Boolean = {
    this.haplotypeString == kp.haplotypeString
  }

  override def compare(otherPath: KmerPath): Int = {
    if (equals(otherPath)) {
      weight.compare(otherPath.weight)
    } else if (weight > otherPath.weight) {
      1
    } else {
      -1
    }
  }
}
