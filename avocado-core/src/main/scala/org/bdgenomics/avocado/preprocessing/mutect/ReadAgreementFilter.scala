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

import htsjdk.samtools.CigarOperator
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rich.ReferenceMappingContext
import org.bdgenomics.avocado.preprocessing.PreprocessingStage
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.adam.rich.RichAlignmentRecord._
import scala.collection.JavaConversions._
import org.bdgenomics.adam.rich.ReferenceMappingContext._

object ReadAgreementFilter extends PreprocessingStage {
  override val stageName: String = "read_agreement_filter"

  private def disagrees(window: Int,
                        rec1: AlignmentRecord, offset1: Int,
                        rec2: AlignmentRecord, offset2: Int): Boolean = {
    (0 until window).exists {
      i =>
        rec1.getSequence.charAt(offset1 + i) != rec2.getSequence.charAt(offset2 + i)
    }
  }

  /**
   * Here's the original language:
   * "If there is an overlapping read pair, and both reads agree the read with the
   * highest quality score is retained otherwise both are discarded."
   *
   * @param rdd
   * @param config
   * @return
   */
  override def apply(rdd: RDD[AlignmentRecord], config: SubnodeConfiguration): RDD[AlignmentRecord] =

    rdd.keyBy(v => ReferenceMappingContext.AlignmentRecordReferenceMapping.getReferenceRegion(v))
      .sortByKey().mapPartitions {

        case itr: Iterator[(ReferenceRegion, AlignmentRecord)] =>
          itr.sliding(2).flatMap {
            case pair: Seq[(ReferenceRegion, AlignmentRecord)] =>
              val (p1, p2) = (pair(0), pair(1))
              val (ref1, rec1) = p1
              val (ref2, rec2) = p2

              /*
            There are four conditions here:
            1. The reads overlap, and
              2. They disagree somewhere along their length -> return neither
              3. They agree along their length -> return the one with higher MAPQ
            4. they don't overlap -> return them both
             */

              // Condition 1
              if (ref1.overlaps(ref2)) {
                val intersect = ref1.intersection(ref2)
                val offset1: Int = intersect.start.toInt - ref1.start.toInt
                val offset2: Int = intersect.start.toInt - ref2.start.toInt

                // Condition 2
                if (disagrees(intersect.length().toInt, rec1, offset1, rec2, offset2)) {
                  Seq()
                } else {
                  // Condition 3
                  if (rec1.getMapq >= rec2.getMapq) {
                    Seq(rec1)
                  } else {
                    Seq(rec2)
                  }
                }

              } else { // Condition 4
                Seq(rec1, rec2)
              }
          }
      }

}
