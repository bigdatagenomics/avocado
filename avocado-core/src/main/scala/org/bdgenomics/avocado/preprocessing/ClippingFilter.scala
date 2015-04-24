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

import htsjdk.samtools.{ Cigar, CigarElement, CigarOperator, TextCigarCodec }
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.formats.avro.AlignmentRecord

object ClippingFilter extends ReadFilter[Double] {
  val stageName: String = "clippingFilter"

  def getThreshold(config: SubnodeConfiguration): Double = {
    // what fraction of bases can be soft clipped?
    val clipped = config.getDouble("fractionSoftClipped", 0.3)
    require(clipped < 1.0 && clipped > 0.0,
      "fractionSoftClipped (%f) must be between 0 and 1, non-inclusive.".format(clipped))
    clipped
  }

  private[preprocessing] def filterRead(r: AlignmentRecord, clipped: Double): Boolean = {
    // what is the sequence length?
    val length = r.getSequence.length.toDouble

    // how many bases are clipped?
    val cigar: List[CigarElement] = TextCigarCodec.decode(r.getCigar)
      .getCigarElements
    val clippedBases = cigar.filter(_.getOperator == CigarOperator.S)
      .map(_.getLength)
      .fold(0)(_ + _)
      .toDouble

    // ratio must be lower than threshold
    clippedBases / length <= clipped
  }
}
