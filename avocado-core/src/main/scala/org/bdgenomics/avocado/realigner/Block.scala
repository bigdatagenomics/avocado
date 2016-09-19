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
package org.bdgenomics.avocado.realigner

import org.bdgenomics.avocado.models.{ Match, ObservationOperator }

/**
 * Sealed trait for describing possibly variant alignment blocks.
 */
private[realigner] sealed trait Block {

  /**
   * Converts this block into alignment operators, possibly applying a function.
   *
   * If this block is not variant, the alignment operator is directly returned.
   * If the block is variant, the function provided is applied to the ref/alt
   * sequences.
   *
   * @param fn A function that takes ref and alt strings and returns alignments.
   * @return Returns the true alignment for this block.
   */
  def fold(fn: (String, String) => Iterable[ObservationOperator]): Iterable[ObservationOperator]
}

/**
 * An alignment block that is an exact sequence match.
 *
 * @param length The length of this sequence match.
 */
private[realigner] case class MatchBlock(length: Int) extends Block {

  /**
   * Directly emits this block as a sequence match.
   *
   * In this fold, the function is not applied.
   *
   * @param fn A function that takes ref and alt strings and returns alignments.
   * @return Returns the true alignment for this block.
   */
  def fold(fn: (String, String) => Iterable[ObservationOperator]): Iterable[ObservationOperator] = {
    Iterable(Match(length))
  }
}

/**
 * An alignment block where the ref and alt sequence do not match.
 *
 * @param ref The reference allele sequence.
 * @param alt The alternate allele sequence.
 */
private[realigner] case class UnknownBlock(ref: String,
                                           alt: String) extends Block {

  assert(ref != alt,
    "Reference (%s) and alternate (%s) sequences match.".format(ref, alt))

  /**
   * Converts this block into alignment operators with a realignment function.
   *
   * @param fn A function that takes ref and alt strings and returns alignments.
   * @return Returns the true alignment for this block.
   */
  def fold(fn: (String, String) => Iterable[ObservationOperator]): Iterable[ObservationOperator] = {
    fn(ref, alt)
  }
}
