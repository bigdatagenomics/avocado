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
package org.bdgenomics.avocado.stats

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.SequenceDictionary
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.formats.avro.{ AlignmentRecord, NucleotideContigFragment }

class AvocadoConfigAndStats(val sc: SparkContext,
                            val debug: Boolean,
                            inputDataset: RDD[AlignmentRecord],
                            val reference: RDD[NucleotideContigFragment]) {

  lazy val coverage = ComputingCoverage.time {
    ScoreCoverage(inputDataset)
  }

  lazy val contigLengths = ReferenceLengths.time {
    GetReferenceContigLengths(reference)
  }

  lazy val referenceSeq = CollectingReference.time {
    reference.collect()
  }

  lazy val sequenceDict = ExtractingSequenceDictionary.time {
    reference.adamGetSequenceDictionary()
  }

  lazy val samplesInDataset = CollectingSamples.time {
    inputDataset.map(_.getRecordGroupSample)
      .distinct()
      .collect()
  }

  lazy val referenceObservations = ExploringReference.time {
    SliceReference(reference)
  }
}
