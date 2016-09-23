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
package org.bdgenomics.avocado.models

import org.scalatest.FunSuite

class ObservationSuite extends FunSuite {

  test("cannot create an observation with empty likelihoods") {
    intercept[AssertionError] {
      Observation(0, 0.0, Array.empty, Array.empty, 1)
    }
  }

  test("cannot create an observation with mismatching likelihood lengths") {
    intercept[AssertionError] {
      Observation(0, 0.0, Array.fill(3)(0.0), Array.fill(2)(0.0), 1)
    }
  }

  test("forward strand must be >= 0") {
    intercept[AssertionError] {
      Observation(-1, 0.0, Array(0.0), Array(0.0), 1)
    }
  }

  test("square map-q must be >= 0") {
    intercept[AssertionError] {
      Observation(0, -1.0, Array(0.0), Array(0.0), 1)
    }
  }

  test("coverage is strictly positive") {
    intercept[AssertionError] {
      Observation(0, 0.0, Array(0.0), Array(0.0), 0)
    }
  }

  test("merge two observations") {
    val obs1 = Observation(0, 0.0, Array(1.0), Array(3.0), 1)
    val obs2 = Observation(1, 1.0, Array(2.0), Array(4.0), 3)

    val newObs = obs1.merge(obs2)

    assert(newObs.forwardStrand === 1)
    assert(newObs.squareMapQ > 0.999 && newObs.squareMapQ < 1.001)
    assert(newObs.alleleLogLikelihoods.size === 1)
    assert(newObs.alleleLogLikelihoods.head > 2.999 && newObs.alleleLogLikelihoods.head < 3.001)
    assert(newObs.otherLogLikelihoods.head > 6.999 && newObs.otherLogLikelihoods.head < 9.001)
    assert(newObs.coverage === 4)
  }

  test("cannot merge two observations that have a different copy number") {
    val obs1 = Observation(0, 0.0, Array(1.0), Array(3.0), 1)
    val obs2 = Observation(1, 1.0, Array(2.0, 2.0), Array(4.0, 4.0), 3)

    intercept[AssertionError] {
      obs1.merge(obs2)
    }
  }
}
