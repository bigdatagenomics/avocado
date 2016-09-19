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

import org.bdgenomics.avocado.models.{
  Deletion,
  Insertion,
  Match,
  ObservationOperator
}
import org.scalatest.FunSuite

class BlockSuite extends FunSuite {

  def foldFn(ref: String, alt: String): Seq[ObservationOperator] = {
    Seq(Deletion(ref), Insertion(alt.length))
  }

  test("folding over a match block returns a match operator") {
    val folded = MatchBlock(5).fold(foldFn).toSeq

    assert(folded.size === 1)
    assert(folded(0) === Match(5))
  }

  test("an unknown block must have mismatching input sequences") {
    intercept[AssertionError] {
      UnknownBlock("AC", "AC")
    }
  }

  test("folding over an unknown block returns our function's result") {
    val folded = UnknownBlock("ACAC", "GT").fold(foldFn).toSeq

    assert(folded.size === 2)
    assert(folded(0) === Deletion("ACAC"))
    assert(folded(1) === Insertion(2))
  }
}
