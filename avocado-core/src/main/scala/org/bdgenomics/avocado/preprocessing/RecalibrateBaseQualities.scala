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
package org.bdgenomics.avocado.preprocessing

import java.io.File
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.SnpTable
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.read.ADAMAlignmentRecordContext._
import org.bdgenomics.formats.avro.AlignmentRecord

object RecalibrateBaseQualities extends PreprocessingStage {

  val stageName = "recalibrateBaseQualities"

  def apply(rdd: RDD[AlignmentRecord], config: SubnodeConfiguration): RDD[AlignmentRecord] = {

    val sc = rdd.sparkContext
    // check for snp table
    val snpTable = if (config.containsKey("snpTable")) {
      sc.broadcast(SnpTable(new File(config.getString("snpTable"))))
    } else {
      sc.broadcast(SnpTable())
    }

    // run bqsr with snp table loaded
    rdd.adamBQSR(snpTable)
  }

}
