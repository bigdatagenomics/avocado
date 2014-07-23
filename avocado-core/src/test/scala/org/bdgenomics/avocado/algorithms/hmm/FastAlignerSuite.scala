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
package org.bdgenomics.avocado.algorithms.hmm

import org.bdgenomics.adam.util.PhredUtils
import org.scalatest.FunSuite
import scala.math.{ abs, log10, pow }

class FastAlignerSuite extends FunSuite {

  def fpEquals(a: Double, b: Double, threshold: Double = 1e-3): Boolean = {
    abs(a - b) <= threshold
  }

  val fa = new FastAligner("ACACTGCACA", 4)

  test("should see hit for exact match") {
    val alignment = fa.alignSequences("ACACTGCACA", "ACAC", "++++")

    assert(alignment.alignmentStateSequence === "MMMMPPPPPP")
    assert(alignment.alignedSequence === "ACAC______")
    assert(alignment.alignedReference === "ACACTGCACA")
    assert(fpEquals(alignment.likelihood, log10(pow(0.9, 4))))
  }

  test("should see match for one error") {
    val alignment = fa.alignSequences("ACACTGCACA", "GCTC", "++++")

    assert(alignment.alignmentStateSequence === "PPPPPMMXMP")
    assert(alignment.alignedSequence === "_____GCTC_")
    assert(alignment.alignedReference === "ACACTGCACA")
    assert(fpEquals(alignment.likelihood, log10(pow(0.9, 3) * 0.1)))
  }

  test("should see exception if we have more than one mismatch") {
    intercept[IllegalArgumentException] {
      val alignment = fa.alignSequences("ACACTGCACA", "GATC", "++++")
    }
  }
}
