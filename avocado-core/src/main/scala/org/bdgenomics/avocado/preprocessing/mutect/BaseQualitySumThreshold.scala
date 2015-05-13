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
package org.bdgenomics.avocado.preprocessing.mutect

import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.avocado.preprocessing.PreprocessingStage
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.adam.rich.RichAlignmentRecord._

object BaseQualitySumThreshold extends PreprocessingStage {
  override val stageName: String = "base_quality_sum_filter"

  override def apply(rdd: RDD[AlignmentRecord], config: SubnodeConfiguration): RDD[AlignmentRecord] =
    rdd.filter {
      case record =>
        record.getMismatchingPositions.toString
          .zip(record.qualityScores).filter(_._1 == 'D').map(_._2).sum <= 100
    }
}
