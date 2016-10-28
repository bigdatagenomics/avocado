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
      Observation(0, 0, 0.0, Array.empty, Array.empty, 1, 0)
    }
  }

  test("cannot create an observation with 1-length likelihoods") {
    intercept[AssertionError] {
      Observation(0, 0, 0.0, Array(0.0), Array(0.0), 1, 0)
    }
  }

  test("cannot create an observation with mismatching likelihood lengths") {
    intercept[AssertionError] {
      Observation(0, 0, 0.0, Array.fill(3)(0.0), Array.fill(2)(0.0), 1, 0)
    }
  }

  test("forward strand must be >= 0") {
    intercept[AssertionError] {
      Observation(-1, 0, 0.0, Array(0.0, 0.0), Array(0.0, 0.0), 1, 0)
    }
    intercept[AssertionError] {
      Observation(0, -1, 0.0, Array(0.0, 0.0), Array(0.0, 0.0), 1, 0)
    }
  }

  test("forward strand cannot exceed coverage") {
    intercept[AssertionError] {
      Observation(2, 0, 0.0, Array(0.0, 0.0), Array(0.0, 0.0), 1, 0)
    }
    intercept[AssertionError] {
      Observation(0, 3, 0.0, Array(0.0, 0.0), Array(0.0, 0.0), 1, 2)
    }
  }

  test("square map-q must be >= 0") {
    intercept[AssertionError] {
      Observation(0, 0, -1.0, Array(0.0, 0.0), Array(0.0, 0.0), 1, 0)
    }
  }

  test("coverage is strictly positive") {
    intercept[AssertionError] {
      Observation(0, 0, 0.0, Array(0.0, 0.0), Array(0.0, 0.0), 0, 0,
        totalCoverage = 0)
    }
    intercept[AssertionError] {
      Observation(0, 0, 0.0, Array(0.0, 0.0), Array(0.0, 0.0), -2, 4)
    }
    intercept[AssertionError] {
      Observation(0, 0, 0.0, Array(0.0, 0.0), Array(0.0, 0.0), 2, -1)
    }
  }

  test("invert an observation") {
    val obs = Observation(0, 1, 0.0, Array(1.0, 2.0), Array(3.0, 4.0), 1, 2,
      totalCoverage = 4)

    val invertedObs = obs.invert
    assert(invertedObs.alleleForwardStrand === 1)
    assert(invertedObs.otherForwardStrand === 0)
    assert(invertedObs.squareMapQ === 0.0)
    assert(invertedObs.copyNumber === 1)
    assert(invertedObs.alleleLogLikelihoods(0) === 3.0)
    assert(invertedObs.alleleLogLikelihoods(1) === 4.0)
    assert(invertedObs.otherLogLikelihoods(0) === 1.0)
    assert(invertedObs.otherLogLikelihoods(1) === 2.0)
    assert(invertedObs.alleleCoverage === 2)
    assert(invertedObs.otherCoverage === 1)
    assert(invertedObs.totalCoverage === 4)
    assert(!invertedObs.isRef)
  }

  test("null an observation") {
    val obs = Observation(0, 1, 0.0, Array(1.0, 2.0), Array(3.0, 4.0), 1, 2,
      totalCoverage = 4)

    val nulledObs = obs.nullOut
    assert(nulledObs.alleleForwardStrand === 0)
    assert(nulledObs.otherForwardStrand === 0)
    assert(nulledObs.squareMapQ === 0.0)
    assert(nulledObs.copyNumber === 1)
    assert(nulledObs.alleleLogLikelihoods(0) === 0.0)
    assert(nulledObs.alleleLogLikelihoods(1) === 0.0)
    assert(nulledObs.otherLogLikelihoods(0) === 1.0)
    assert(nulledObs.otherLogLikelihoods(1) === 2.0)
    assert(nulledObs.alleleCoverage === 0)
    assert(nulledObs.otherCoverage === 0)
    assert(nulledObs.totalCoverage === 4)
    assert(!nulledObs.isRef)
  }

  test("merge two observations") {
    val obs1 = Observation(0, 1, 0.0, Array(1.0, 2.0), Array(3.0, 4.0), 1, 2,
      totalCoverage = 3)
    val obs2 = Observation(1, 2, 1.0, Array(2.0, 4.0), Array(4.0, 6.0), 3, 3,
      totalCoverage = 7)

    val newObs = obs1.merge(obs2)

    assert(newObs.alleleForwardStrand === 1)
    assert(newObs.otherForwardStrand === 3)
    assert(newObs.squareMapQ > 0.999 && newObs.squareMapQ < 1.001)
    assert(newObs.copyNumber === 1)
    assert(newObs.alleleLogLikelihoods(0) > 2.999 && newObs.alleleLogLikelihoods(0) < 3.001)
    assert(newObs.alleleLogLikelihoods(1) > 5.999 && newObs.alleleLogLikelihoods(1) < 6.001)
    assert(newObs.otherLogLikelihoods(0) > 6.999 && newObs.otherLogLikelihoods(0) < 7.001)
    assert(newObs.otherLogLikelihoods(1) > 9.999 && newObs.otherLogLikelihoods(1) < 10.001)
    assert(newObs.alleleCoverage === 4)
    assert(newObs.otherCoverage === 5)
    assert(newObs.coverage === 9)
    assert(newObs.totalCoverage === 10)
    assert(newObs.isRef)
  }

  test("merge two alt observations") {
    val obs1 = Observation(0, 1, 0.0, Array(1.0, 2.0), Array(3.0, 4.0), 1, 2,
      totalCoverage = 3, isRef = false)
    val obs2 = Observation(1, 2, 1.0, Array(2.0, 4.0), Array(4.0, 6.0), 3, 3,
      totalCoverage = 7, isRef = false)

    val newObs = obs1.merge(obs2)

    assert(newObs.alleleForwardStrand === 1)
    assert(newObs.otherForwardStrand === 3)
    assert(newObs.squareMapQ > 0.999 && newObs.squareMapQ < 1.001)
    assert(newObs.copyNumber === 1)
    assert(newObs.alleleLogLikelihoods(0) > 2.999 && newObs.alleleLogLikelihoods(0) < 3.001)
    assert(newObs.alleleLogLikelihoods(1) > 5.999 && newObs.alleleLogLikelihoods(1) < 6.001)
    assert(newObs.otherLogLikelihoods(0) > 6.999 && newObs.otherLogLikelihoods(0) < 7.001)
    assert(newObs.otherLogLikelihoods(1) > 9.999 && newObs.otherLogLikelihoods(1) < 10.001)
    assert(newObs.alleleCoverage === 4)
    assert(newObs.otherCoverage === 5)
    assert(newObs.coverage === 9)
    assert(newObs.totalCoverage === 10)
    assert(!newObs.isRef)
  }

  test("cannot merge observations that disagree about whether the allele is an alt") {
    val obs1 = Observation(0, 1, 0.0, Array(1.0, 2.0), Array(3.0, 4.0), 1, 2,
      totalCoverage = 3)
    val obs2 = Observation(1, 2, 1.0, Array(2.0, 4.0), Array(4.0, 6.0), 3, 3,
      totalCoverage = 7, isRef = false)

    intercept[AssertionError] {
      obs1.merge(obs2)
    }
  }

  test("cannot merge two observations that have a different copy number") {
    val obs1 = Observation(0, 1, 0.0, Array(1.0, 1.0), Array(3.0, 3.0), 1, 2)
    val obs2 = Observation(1, 1, 1.0, Array(2.0, 2.0, 2.0), Array(4.0, 4.0, 4.0), 3, 2)

    intercept[AssertionError] {
      obs1.merge(obs2)
    }
  }
}
