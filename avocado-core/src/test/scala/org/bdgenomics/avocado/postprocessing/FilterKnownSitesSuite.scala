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
package org.bdgenomics.avocado.postprocessing

import org.bdgenomics.adam.models.SnpTable
import org.bdgenomics.adam.rich.RichVariant._
import org.bdgenomics.formats.avro.{ Contig, Variant }
import org.scalatest.FunSuite

class FilterKnownSitesSuite extends FunSuite {

  val snpTable = new SnpTable(Map("chr1" -> Set(10000L)))
  val keepFilter = new FilterKnownSites(None, snpTable, true)
  val discardFilter = new FilterKnownSites(None, snpTable, false)
  val ctg = Contig.newBuilder()
    .setContigName("chr1")
    .build()

  test("process a known site correctly") {
    val v = Variant.newBuilder()
      .setContig(ctg)
      .setStart(10000L)
      .build()
    assert(keepFilter.filterFn(v))
    assert(!discardFilter.filterFn(v))
  }

  test("process an unknown site correctly") {
    val v = Variant.newBuilder()
      .setContig(ctg)
      .setStart(10001L)
      .build()
    assert(!keepFilter.filterFn(v))
    assert(discardFilter.filterFn(v))
  }
}
