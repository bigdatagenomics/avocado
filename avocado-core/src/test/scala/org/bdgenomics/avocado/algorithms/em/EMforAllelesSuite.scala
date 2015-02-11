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
package org.bdgenomics.avocado.algorithms.em

import org.bdgenomics.avocado.algorithms.math.MathTestUtils
import org.scalatest.FunSuite
import scala.math.{ abs, exp, log }

class EMForAllelesSuite extends FunSuite {

  test("cannot run EM without specifying an iteration limit or target tolerance") {
    intercept[AssertionError] {
      EMForAlleles.emForMAF(Array(Array[Double]()),
        1.0 - 1e-3)
    }
  }

  test("run EM on single sample, definite ref") {
    val psi = EMForAlleles.emForMAF(Array(Array(1e-12, 1e-6, 1.0 - 1e-6 - 1e-12).map(log(_))),
      log(1.0 - 1e-3),
      maxIterations = Some(10))

    MathTestUtils.assertAlmostEqual(exp(psi), 1.0)
  }

  test("run EM on three samples, mix of hom ref, het, hom alt") {
    val psi = EMForAlleles.emForMAF(Array(Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.0001, 0.1, 0.0001).map(log(_)),
      Array(0.1, 0.0001, 0.0000001).map(log(_))),
      log(1.0 - 1e-3),
      maxIterations = Some(10))

    MathTestUtils.assertAlmostEqual(exp(psi), 0.5)
  }

  test("run EM on five samples, one hom alt, all others hom ref") {
    val psi = EMForAlleles.emForMAF(Array(Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.1, 0.0001, 0.0000001).map(log(_))),
      log(1.0 - 1e-3),
      maxIterations = Some(10))

    MathTestUtils.assertAlmostEqual(exp(psi), 0.8)
  }

  test("run EM on five samples, with varying ploidy, M = 10, G = 7") {
    val psi = EMForAlleles.emForMAF(Array(Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.0000001, 0.0001, 0.1).map(log(_)),
      Array(0.1, 0.0001).map(log(_)),
      Array(0.0001, 0.1, 0.0001, 0.0000001).map(log(_))),
      log(1.0 - 1e-3),
      maxIterations = Some(10))

    MathTestUtils.assertAlmostEqual(exp(psi), 0.7)
  }
}
