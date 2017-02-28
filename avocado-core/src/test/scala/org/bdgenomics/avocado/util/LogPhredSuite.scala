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
package org.bdgenomics.avocado.util

import org.scalatest.FunSuite
import scala.math.log

class LogPhredSuite extends FunSuite {

  test("convert log error probabilities to phred scores") {
    val phred10 = LogPhred.logErrorToPhred(log(0.1))
    assert(phred10 > 9.999 && phred10 < 10.001)
    val phred50 = LogPhred.logErrorToPhred(log(0.00001))
    assert(phred50 > 49.999 && phred50 < 50.001)
  }
}
