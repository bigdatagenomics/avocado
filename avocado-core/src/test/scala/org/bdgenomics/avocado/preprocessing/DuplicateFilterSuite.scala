/*
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.avocado.preprocessing

import org.bdgenomics.formats.avro.AlignmentRecord
import org.scalatest.FunSuite

class DuplicateFilterSuite extends FunSuite {

  test("should filter out a read that is marked as a duplicate") {
    val read = AlignmentRecord.newBuilder()
      .setDuplicateRead(true)
      .build()

    assert(!DuplicateFilter.filterRead(read, true))
  }

  test("should keep a unique read") {
    val read = AlignmentRecord.newBuilder()
      .setDuplicateRead(false)
      .build()

    assert(DuplicateFilter.filterRead(read, true))
  }

  test("should keep an unmarked read") {
    val read = AlignmentRecord.newBuilder()
      .build()

    assert(DuplicateFilter.filterRead(read, true))
  }
}
