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

import org.bdgenomics.adam.models.ReferencePosition

/**
 * A kmer has a prefix of length k - 1 and a unit length suffix.
 *
 * @param kmerSeq sequence that will be split into k-1 len prefix and single char suffix
 * @param weight
 */
private[debrujin] case class Kmer(kmerSeq: String,
                                  refPos: Option[ReferencePosition] = None,
                                  var phred: List[Int] = List(),
                                  var mapq: List[Int] = List(),
                                  var readId: List[Long] = List(),
                                  var isNegativeStrand: List[Boolean] = List(),
                                  var predecessors: List[Kmer] = List(),
                                  var successors: List[Kmer] = List()) {

  assert(phred.length == mapq.length)

  def multiplicity = phred.length
  val prefix: String = kmerSeq.dropRight(1)
  val suffix: Char = kmerSeq.last
  val isReference = refPos.isDefined

  def nextPrefix: String = {
    prefix.drop(1) + suffix
  }

  override def toString: String = {
    prefix + "[" + suffix + "]"
  }

  def toDetailedString: String = {
    kmerSeq + " * " + multiplicity + ", " + refPos.fold("unmapped")("@ " + _) + "\n" +
      "qual: " + phred.map(_.toString).fold("")(_ + ", " + _) + "\n" +
      "mapq: " + mapq.map(_.toString).fold("")(_ + ", " + _) + "\n" +
      "readId: " + readId.map(_.toString).fold("")(_ + ", " + _) + "\n" +
      "pre: " + predecessors.map(_.toString).fold("")(_ + ", " + _) + "\n" +
      "post: " + successors.map(_.toString).fold("")(_ + ", " + _)
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

  def removeLinks(kmer: Kmer) {
    predecessors = predecessors.filter(_.kmerSeq != kmer.kmerSeq)
    successors = successors.filter(_.kmerSeq != kmer.kmerSeq)
  }
}
