/*
 * Copyright (c) 2014. Regents of the University of California.
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

package org.bdgenomics.avocado.algorithms.hmm

import org.scalatest.FunSuite
import scala.math.{ log, log10, pow }

class HaplotypePairArgmaxSuite extends FunSuite {

  test("apply softmax function to log scaled numbers") {
    val aIn = Array(log(0.1), log(0.3))
    val aOut = HaplotypePairArgmax.softmax(aIn)

    assert(aOut(0) > 0.24 && aOut(0) < 0.26)
    assert(aOut(1) > 0.74 && aOut(1) < 0.76)
  }

  test("apply softmax function to integers") {
    val a = HaplotypePairArgmax.softmax((0 to 5).toArray.map(_.toDouble))
    val s = a.sum

    assert(s > 0.99 && s < 1.01)
  }

  test("throw assert for invalid input sizes passed to softmax") {
    intercept[AssertionError] {
      HaplotypePairArgmax.softmax(Array[Double]())
    }
  }

  test("apply ratio") {
    val aIn = Array(1.0, 3.0)
    val aOut = HaplotypePairArgmax.ratio(aIn)

    assert(aOut(0) > 0.24 && aOut(0) < 0.26)
    assert(aOut(1) > 0.74 && aOut(1) < 0.76)
  }

  test("throw assert for invalid formats passed to ratio") {
    intercept[AssertionError] {
      HaplotypePairArgmax.ratio(Array[Double]())
    }

    intercept[AssertionError] {
      HaplotypePairArgmax.ratio(Array[Double](0.1, -0.1))
    }
  }

  test("throw assert for invalid input sizes passed to EM algorithm") {
    intercept[AssertionError] {
      HaplotypePairArgmax(Array(1L, 1L), Array(Array(1.0), Array(1.0), Array(1.0)), 0, 0.0, 0.0)
    }

    intercept[AssertionError] {
      HaplotypePairArgmax(Array[Long](), Array[Array[Double]](), 0, 0.0, 0.0)
    }

    intercept[AssertionError] {
      HaplotypePairArgmax(Array(1L, 1L), Array(Array(1.0), Array(1.0, 1.0)), 0, 0.0, 0.0)
    }
  }

  test("throw assert for invalid input sizes passed to transpose") {
    intercept[AssertionError] {
      HaplotypePairArgmax.transpose(Array[Array[Double]]())
    }

    intercept[AssertionError] {
      HaplotypePairArgmax.transpose(Array(Array[Double]()))
    }

    intercept[AssertionError] {
      HaplotypePairArgmax.transpose(Array(Array(0.0), Array(1.0, 0.0)))
    }
  }

  test("transpose single row/column matrices") {
    val row = (0 to 5).toArray.map(i => Array(i.toDouble))
    val tcol = HaplotypePairArgmax.transpose(row)

    assert(tcol.length === 1)
    assert(tcol.forall(_.length == 6))
    assert((0 to 5).forall(i => {
      tcol(0)(i) == i.toDouble
      tcol(0)(i) == row(i)(0)
    }))

    val col = Array((0 to 5).toArray.map(_.toDouble))
    val trow = HaplotypePairArgmax.transpose(col)

    assert(trow.length === 6)
    assert(trow.forall(_.length == 1))
    assert((0 to 5).forall(i => {
      trow(i)(0) == i.toDouble
      trow(i)(0) == col(0)(i)
    }))
  }

  test("compute factorials") {
    assert(1 === pow(10.0, HaplotypePairArgmax.log10Factorial(0)).toInt)
    assert(1 === pow(10.0, HaplotypePairArgmax.log10Factorial(1)).toInt)
    assert(6 === pow(10.0, HaplotypePairArgmax.log10Factorial(3)).toInt)
    assert(720 === pow(10.0, HaplotypePairArgmax.log10Factorial(6)).round.toInt)
    assert((1L to 20L).reduce(_ * _) * 0.99999 < pow(10.0, HaplotypePairArgmax.log10Factorial(20)))
    assert((1L to 20L).reduce(_ * _) * 1.00001 > pow(10.0, HaplotypePairArgmax.log10Factorial(20)))
  }

  test("compute multinomial for one trial") {
    val pTrial = HaplotypePairArgmax.multinomial(Array(0.4, 0.6), Array(1))

    // multinomial returns log10 scored values
    assert(pTrial > log10(0.59))
    assert(pTrial < log10(0.61))
  }

  test("compute multinomial for two trials") {
    val pTrial = HaplotypePairArgmax.multinomial(Array(0.4, 0.35, 0.25), Array(0, 2))

    // multinomial returns log10 scored values
    assert(pTrial > log10(0.19))
    assert(pTrial < log10(0.21))
  }

  test("check asserts for multinomial") {
    intercept[AssertionError] {
      HaplotypePairArgmax.multinomial(Array(0.0, 0.0), Array(0))
    }

    intercept[AssertionError] {
      HaplotypePairArgmax.multinomial(Array(1.0, 1.0), Array(0))
    }

    intercept[AssertionError] {
      HaplotypePairArgmax.multinomial(Array(0.0, 1.0), Array())
    }

    intercept[AssertionError] {
      HaplotypePairArgmax.multinomial(Array(0.0, 1.0), Array(5))
    }
  }

  test("pick given threshold") {
    assert(HaplotypePairArgmax.pick(0.4, Array(0.45, 0.55)) === 0)

    val picks = HaplotypePairArgmax.pickReads(0.4, Array(Array(0.45, 0.55), Array(0.35, 0.65)))
    assert(picks(0) === 0)
    assert(picks(1) === 1)
  }

  test("choose best read mappings for simple case") {
    // equal length haplotypes
    val lengths = Array(10L, 10L)
    // all reads map unambiguously
    val mapping = Array(Array(-0.0, -0.0, -10.0, -10.0),
      Array(-10.0, -10.0, -0.0, -0.0))

    // generate assignments
    val (l, assignments) = HaplotypePairArgmax.run(lengths, mapping, 10, 0.05, 0.5)

    // check obvious assignments
    assert(assignments(0) === 0)
    assert(assignments(1) === 0)
    assert(assignments(2) === 1)
    assert(assignments(3) === 1)

    assert(l > -0.426 && l < -0.425)
  }

  test("choose best read mappings for more complex case") {
    // diffrent length haplotypes
    val lengths = Array(6L, 4L)
    // all but one read are mapped unambiguously
    val mapping = Array(Array(-0.0, -0.0, -10.0, -10.0, -4.0),
      Array(-10.0, -10.0, -0.0, -0.0, -3.9))

    // generate assignments
    val (l, assignments) = HaplotypePairArgmax.run(lengths, mapping, 10, 0.05, 0.5)

    // check obvious assignments
    assert(assignments(0) === 0)
    assert(assignments(1) === 0)
    assert(assignments(2) === 1)
    assert(assignments(3) === 1)
    assert(assignments(4) === 0)
  }
}
