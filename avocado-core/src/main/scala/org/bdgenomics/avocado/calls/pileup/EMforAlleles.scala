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
package org.bdgenomics.avocado.calls.pileup

import scala.math.pow

object EMForAlleles {

  /**
   * Main MAF EM function
   *  IN: Phi - an initial MAF vector of length number of SNPs
   *      GL - Array of arrays of likelihood triples P( D | g )
   *          (note these are NOT multiplied by P(g | phi)! )
   *  OUT: Phi - ML estimate of MAF's across SNPs
   *
   *  Note: GL is currently an array of (numSnps) arrays of length (numInds),
   */
  def emForMAF(Phi: Array[Double], GL: Array[Array[(Double, Double, Double)]]): Array[Double] = {
    var eps = 1.0
    val tol = 0.0001
    val L = Phi.length
    var phi_updates = Phi
    while (eps > tol) {
      var Phi_next = Array.fill(L) { 0.0 }
      Phi_next.indices.foreach(i => {
        GL(i).foreach(l => {
          Phi_next(i) += (1.0 / (2.0 * GL(i).length)) * ((1.0 * l._2 * 2.0 * phi_updates(i) * (1 - phi_updates(i)) +
            2.0 * l._3 * pow(phi_updates(i), 2.0)) / (l._1 * pow(1.0 - phi_updates(i), 2.0) +
              l._2 * 2.0 * phi_updates(i) * (1.0 - phi_updates(i)) +
              l._3 * pow(phi_updates(i), 2.0)))
        })
      })
      var eps = 0.0
      phi_updates.indices.foreach(i => eps += pow(phi_updates(i) - Phi_next(i), 2.0))
      phi_updates = Phi_next
    }
    return phi_updates
  }

  /**
   * Helper function to compute Y iteratively
   * For each site, executes the recursion in 4.2.3. Y(i) is Ynk vector for site i
   */
  def compY(GL: Array[Array[(Double, Double, Double)]]): Array[Array[Double]] = {
    val L = GL.length
    val GLt = GL.transpose
    val n = GLt.length
    val M = 2 * n
    var Y = Array.ofDim[Double](L, n + 1, M + 1)
    // NOTE: this ordering may be suboptimal?
    for (i <- 0 until L) {
      for (k <- 0 to M) {
        for (j <- 0 to n) { // 0 = 0 people not first person
          if (j == 0) {
            Y(i)(j)(k) = 1.0
          } else if (k == 0) {
            Y(i)(j)(k) = (1.0 / (2.0 * j * (2.0 * j - 1.0))) * ((2.0 * j - k) * (2.0 * j - k - 1.0) * Y(i)(j - 1)(k) * GL(i)(j)._1)
          } else if (k == 1) {
            Y(i)(j)(k) = (1.0 / (2.0 * j * (2.0 * j - 1.0))) * ((2.0 * j - k) * (2.0 * j - k - 1.0) * Y(i)(j - 1)(k) * GL(i)(j)._1 +
              2.0 * k * (2.0 * j - k) * Y(i)(j - 1)(k - 1) * GL(i)(j)._2)
          } else {
            Y(i)(j)(k) = (1.0 / (2.0 * j * (2.0 * j - 1.0))) * ((2.0 * j - k) * (2.0 * j - k - 1.0) * Y(i)(j - 1)(k) * GL(i)(j)._1 +
              2.0 * k * (2.0 * j - k) * Y(i)(j - 1)(k - 1) * GL(i)(j)._2 + k * (k - 1.0) *
              Y(i)(j - 1)(k - 2) * GL(i)(j)._2)
          }
        }
      }
    }

    var Yr = Array.ofDim[Double](L, M)
    for (l <- 0 until L) Yr(l) = Y(l)(n)
    return Yr
  }

  /**
   * Main AFS EM function
   *   IN: Phi - an initial MAF vector of length number of SNPs
   *       GL - Array of arrays of likelihood triples P( D | g )
   *           (note these are NOT multiplied by P(g | phi)! )
   *   OUT: Phi - ML estimate of MAF's across SNPs
   *   Note: GL is currently an array of (numSnps) arrays of length (numInds), which is transposed
   */
  def emForAFS(Phik: Array[Double], GL: Array[Array[(Double, Double, Double)]]): Array[Double] = {
    val GLt = GL.transpose
    val tol = 0.0001
    val L = GL.length
    val M = Phik.length
    var eps = 1.0
    var Y = compY(GL)
    var phik_updates = Phik
    while (eps > tol) {
      var sums = Array.fill(L) { 0.0 }
      sums.indices.foreach(a => phik_updates.indices.foreach(p => sums(a) += phik_updates(p) * Y(a)(p)))
      val Phik_next = Array.fill(M) { 0.0 }
      Phik_next.indices.foreach(i => Y.foreach(y => Phik_next(i) += (1.0 / L) * phik_updates(i) * y(i) / sums(i)))
      eps = 0.0
      phik_updates.indices.foreach(i => eps += pow(phik_updates(i) - Phik_next(i), 2.0))
      phik_updates = Phik_next
    }
    phik_updates
  }
}
