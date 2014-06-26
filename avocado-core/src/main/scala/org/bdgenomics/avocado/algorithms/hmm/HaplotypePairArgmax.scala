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

import scala.math.{ exp, pow, log10 }
import scala.util.Random

private[hmm] object HaplotypePairArgmax {

  /**
   *
   * @param haplotypeLengths An array containing the lengths of each haplotype.
   * @param mappingLikelihoods Predetermined likelihoods of each read mapping to each haplotype. For _n_
   * haplotypes and _k_ reads, matrix should be n x k.
   * @return Returns a tuple containing the haplotype pair likelihood, as well as haplotype assignments per read.
   */
  def apply(haplotypeLengths: Array[Long],
            mappingLikelihoods: Array[Array[Double]],
            maxIters: Int,
            epsilon: Double,
            stepSize: Double): (Double, Array[Int]) = {

    // check haplotype lengths and matrix formulations
    assert(haplotypeLengths.length == mappingLikelihoods.length, "Sizes of input arrays are inconsistent.")
    assert(haplotypeLengths.length > 1, "Input matrices must contain at least 2 haplotypes.")
    assert(mappingLikelihoods.forall(_.length == mappingLikelihoods.head.length),
      "Mapping likelihood matrix rows must have consistent lengths.")

    // run
    run(haplotypeLengths, mappingLikelihoods, maxIters, epsilon, stepSize)
  }

  /**
   * Computes the softmax of an array of n elements. Softmax is defined as:
   *
   * \sigma(v)_z = \frac{\exp{v_z}}{\sum_i \exp{v_i}}
   *
   * And is used to map a vector of real values to a vector that sums to 1, for use in a categorical distribution.
   *
   * @param v Vector to apply softmax to.
   * @return Vector that sums to 1.
   */
  def softmax(v: Array[Double]): Array[Double] = {

    // softmax must have at least 1 element
    assert(v.length > 0, "Cannot run softmax on empty array.")

    val expArray = v.map(exp)
    val norm = expArray.sum

    expArray.map(_ / norm)
  }

  /**
   * Takes the ratio of all elements in a matrix, to the size of the whole.
   *
   * @param v Vector to compute ratio over. All elements must be greater than or equal to 0.
   * @return Vector that sums to 1 containing ratios.
   */
  def ratio(v: Array[Double]): Array[Double] = {

    // ratio must have at least 1 element and all elements must be >= 0
    assert(v.length > 0, "Cannot run ratio on empty array.")
    assert(v.forall(_ >= 0.0), "All elements must be greater than or equal to 0.")

    val t = v.sum

    v.map(_ / t)
  }

  /**
   * Transposes an _n_ * _m_ matrix to an _m_ * _n_ matrix.
   *
   * @param matrix Matrix to transpose.
   * @return Transposed matrix.
   */
  def transpose(matrix: Array[Array[Double]]): Array[Array[Double]] = {

    // check matrix sizes
    assert(matrix.forall(_.length == matrix.head.length), "Matrix rows must all be the same length")
    assert(matrix.length > 0, "Matrix must contain elements.")
    assert(matrix.head.length > 0, "Matrix must contain elements.")

    // allocate matrix rows
    val tMatrix = new Array[Array[Double]](matrix.head.length)

    (0 until matrix.head.length).foreach(i => {

      // allocate new matrix column
      tMatrix(i) = new Array[Double](matrix.length)

      (0 until matrix.length).foreach(j => {
        // copy element
        tMatrix(i)(j) = matrix(j)(i)
      })
    })

    tMatrix
  }

  /**
   * Computes the log10 factorial of a given number. This takes an int, but returns a double.
   */
  def log10Factorial(i: Int): Double = {
    (1 to i).map(n => log10(n.toDouble)).sum
  }

  /**
   * Computes the log10 probability of a multinomial random, given _n_ trials. E.g.,
   *
   *       p(n | theta) =       \frac{n!}{\prod x_i!}   \prod p_i^{x_i}
   * log10 p(n | theta) = log10 \frac{n!}{\prod x_i!} + \sum  x_i log10 p_i
   *
   * @param probabilities Success probabilities per class.
   * @param trial Observed values from _n_ trials.
   * @return The probability of observing the _n_ trials seen before.
   */
  def multinomial(probabilities: Array[Double], trials: Array[Int]): Double = {
    val psum = probabilities.sum
    assert(psum > 0.99999 && psum < 1.00001, "Sum of probabilities must be 1.")
    assert(trials.length > 0, "Must have at least 1 trial.")
    assert(trials.forall(t => t >= 0 && t < probabilities.length), "Trial must be in range.")

    // get occurrence count
    val count = (0 until probabilities.length).map(c => (c, trials.count(_ == c))).filter(p => p._2 != 0)

    // compute probability
    val nFact = log10Factorial(trials.length)
    val per = count.map(p => {
      // get class and count
      val (cl, count) = p

      log10(probabilities(cl)) * count - log10Factorial(count)
    })

    nFact + per.reduce(_ + _)
  }

  /**
   * Picks the class of a read given a likelihood threshold.
   */
  def pick(threshold: Double, likelihoods: Array[Double], r: Random = new Random()): Int = {
    assert(likelihoods.length == 2)

    if (likelihoods(0) == threshold) {
      r.nextInt(2) // needed to break ties, else have wierd scoring for homozygous cases
    } else if (likelihoods(0) > threshold) {
      0
    } else {
      1
    }
  }

  /**
   * Evaluates the log likelihood of a class mapping based on a threshold value.
   *
   * @param threshold Ratio based threshold to apply.
   * @param pReadsPerHaplotype The probability of a single read mapping to a given haplotype.
   * @param rMappingToHaplotype Softmax ratio of mapping likelihoods for all reads.
   * @param lMappingToHaplotype Log scale mapping likelihoods for reads to haplotypes.
   * @return Returns a tuple containing haplotype assignments and log likelihood.
   */
  def evaluate(threshold: Double,
               pReadsPerHaplotype: Array[Double],
               rMappingToHaplotype: Array[Array[Double]],
               lMappingToHaplotype: Array[Array[Double]],
               r: Random = new Random()): (Array[Int], Double) = {

    // get classes
    val classes = pickReads(threshold, rMappingToHaplotype, r)

    // compute likelihoods
    val classLikelihoods: Array[Double] = lMappingToHaplotype.zip(classes).map((v: (Array[Double], Int)) => v._1(v._2))

    (classes, multinomial(pReadsPerHaplotype, classes) + classLikelihoods.reduce(_ + _))
  }

  /**
   * Picks the classes for a set of reads.
   */
  def pickReads(threshold: Double, likelihoods: Array[Array[Double]], r: Random = new Random()): Array[Int] = {
    likelihoods.map(pick(threshold, _, r))
  }

  def run(haplotypeLengths: Array[Long],
          mappingLikelihoods: Array[Array[Double]],
          maxIters: Int,
          epsilon: Double,
          stepSize: Double): (Double, Array[Int]) = {

    val r = new Random

    // we use the haplotype lengths as the mixing proportions
    val proportions = ratio(haplotypeLengths.map(_.toDouble))

    // transpose the likelihoods and apply a softmax
    val transposedMappingLikelihoods = transpose(mappingLikelihoods)
    val normalizedMappingLikelihoods = transposedMappingLikelihoods.map(softmax)

    // initialize
    var step = 0.5
    var point = 0.5
    var iter = 0
    var (classes, l) = evaluate(point,
      proportions,
      normalizedMappingLikelihoods,
      transposedMappingLikelihoods,
      r)
    var improvement = epsilon + 1.0

    // loop and try to find the best location
    do {
      step = stepSize * step

      // try higher and lower
      val (higherClasses, hL) = evaluate(point + step,
        proportions,
        normalizedMappingLikelihoods,
        transposedMappingLikelihoods)
      val (lowerClasses, lL) = evaluate(point - step,
        proportions,
        normalizedMappingLikelihoods,
        transposedMappingLikelihoods)

      // accept/reject step?
      if (hL > l && hL > lL) { // accept higher
        classes = higherClasses
        improvement = hL - l
        point = point + step
        l = hL
      } else if (lL > l) { // accept lower
        classes = lowerClasses
        improvement = lL - l
        point = point - step
        l = lL
      } // reject

      iter += 1
    } while (iter < maxIters && improvement > epsilon && point > 0.0 && point < 1.0)

    (l, classes)
  }

}
