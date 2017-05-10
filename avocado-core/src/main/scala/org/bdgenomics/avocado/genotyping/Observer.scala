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
package org.bdgenomics.avocado.genotyping

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.util.PhredUtils
import org.bdgenomics.avocado.models.{
  Clipped,
  Deletion,
  Insertion,
  Match,
  Observation,
  ObservationOperator
}
import org.bdgenomics.formats.avro.AlignmentRecord
import scala.math.log

/**
 * Produces raw allelic likelihoods. This is how you get [[Observation]]s!
 */
private[genotyping] object Observer extends Serializable {

  /**
   * From a single read, emits likelihood observations.
   *
   * Emits a likelihood for each allele seen in this read.
   *
   * @param read Read to observe.
   * @param ploidy Sample ploidy to assume.
   * @return Returns an iterable collection of observations, keyed by (site,
   *   allele, sample ID).
   */
  def observeRead(read: AlignmentRecord,
                  ploidy: Int): Iterable[((ReferenceRegion, String, String), SummarizedObservation)] = {

    // extract cigar
    val alignment = ObservationOperator.extractAlignmentOperators(read)

    // for convenience, get the sample name, mapping quality, sequence,
    // qualities, and the contig name
    val sampleId = read.getRecordGroupSample
    val contigName = read.getContigName
    val mapQ = read.getMapq
    val readSequence = read.getSequence
    val readQualities = read.getQual
    val forwardStrand = !read.getReadNegativeStrand

    // map over the alignment operators and generate allelic observations
    var readIdx = 0
    var pos = read.getStart
    alignment.flatMap(op => op match {
      case Clipped(length, isSoft) => {
        if (isSoft) {
          readIdx += length
        }
        Iterable.empty
      }
      case Match(length, optRef) => {
        (0 until length).map(idx => {

          // the key is the (site, allele, sampleId)
          val key = (ReferenceRegion(contigName, pos, pos + 1),
            readSequence(readIdx).toString,
            sampleId)

          // build the observation
          val obs = SummarizedObservation(optRef.isEmpty,
            forwardStrand,
            Some(readQualities(readIdx).toInt - 33),
            mapQ)

          // increment the indices
          readIdx += 1
          pos += 1

          (key, obs)
        })
      }
      case Insertion(length) => {

        // store old index and update
        val oldReadIdx = readIdx
        readIdx += length

        // get bases and quals
        val bases = readSequence.substring(oldReadIdx, readIdx)
        val qual = readQualities.substring(oldReadIdx, readIdx)
          .map(_.toInt - 33)
          .sum / length

        // the key is the (site, allele, sampleId)
        // insertions associate to the site to their left, hence the -1
        val key = (ReferenceRegion(contigName, pos - 1, pos),
          bases,
          sampleId)

        // build the observation
        val obs = SummarizedObservation(false,
          forwardStrand,
          Some(qual),
          mapQ)

        Iterable((key, obs))
      }
      case d: Deletion => {

        // copy and update pos
        val oldPos = pos
        pos += d.size

        // the key is the (site, allele, sampleId)
        // deletions have an empty string for the allele
        val key = (ReferenceRegion(contigName, oldPos, pos),
          "",
          sampleId)

        // build the observation
        val obs = SummarizedObservation(false,
          forwardStrand,
          None,
          mapQ)

        Iterable((key, obs))
      }
    })
  }

  /**
   * Compute the allelic likelihoods for a sample, given a fixed copy number.
   *
   * @param copyNumber The assumed copy number at the site.
   * @param mapSuccessProb The probability that this read was mapped to the
   *   correct site.
   * @param baseQuality The optional base quality (in Phred) for this site.
   * @return Returns a tuple of the (allele, non-allele) likelihoods.
   */
  private[genotyping] def likelihoods(copyNumber: Int,
                                      mapSuccessProb: Double,
                                      baseQuality: Option[Int]): (Array[Double], Array[Double]) = {

    // build allele/other likelihood arrays
    val alleleArray = new Array[Double](copyNumber + 1)
    val otherArray = new Array[Double](copyNumber + 1)

    // go base quality to phred
    val baseSuccessProb = baseQuality.fold(1.0)(
      q => PhredUtils.phredToSuccessProbability(q))

    // epsilon and 1 - epsilon
    val oneMinusEpsilon = mapSuccessProb * baseSuccessProb
    val epsilon = 1.0 - oneMinusEpsilon

    // populate arrays with likelihoods
    var g = 0
    while (g <= copyNumber) {
      val mMinusG = copyNumber - g
      alleleArray(g) = log(mMinusG * epsilon + g * oneMinusEpsilon)
      otherArray(g) = log(mMinusG * oneMinusEpsilon + g * epsilon)
      g += 1
    }

    (alleleArray, otherArray)
  }
}
