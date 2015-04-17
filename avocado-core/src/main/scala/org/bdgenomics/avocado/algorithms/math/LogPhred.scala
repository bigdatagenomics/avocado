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
package org.bdgenomics.avocado.algorithms.math

import scala.math.{ expm1, log => mathLog }

object LogToPhred {

  private val LOG10 = mathLog(10.0)

  /**
   * Conversion between log _success_ probabilities and Phred scores.
   *
   * Q = -10 log_{10} (1 - p)
   *
   * If l = ln p, then:
   *
   * l_e = ln 1 + ln (1 - exp(l_e / 1))
   *     = ln(1 - exp(l_e))
   * Q = -10 l_e / ln 10
   *
   * @note Just remember kids: if you want to live a happy, fulfilling life, don't use Phred!
   * @note We don't have a conversion the other way around because it doesn't really make sense.
   *
   * @param l Log probability.
   * @return Returns the Phred score corresponding to that log probability.
   */
  def log2phred(l: Double): Double = {
    (-10.0 * mathLog(-expm1(l)) / LOG10)
  }
}
