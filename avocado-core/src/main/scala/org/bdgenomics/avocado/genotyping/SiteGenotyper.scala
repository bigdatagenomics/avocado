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

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.{ InstrumentedRDD, RDD }
import org.bdgenomics.adam.models._
import org.bdgenomics.adam.rdd.GenomicPositionPartitioner
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.formats.avro._
import scala.math.max

trait SiteGenotyper extends Genotyper {

  val contigLengths: Map[String, Long]

  final def genotype(observations: RDD[Observation]): RDD[VariantContext] = {
    val sortedRdd = observations.keyBy(_.pos)
      .repartitionAndSortWithinPartitions(GenomicPositionPartitioner(observations.partitions.size,
        contigLengths))
    val instrumentedRdd = {
      import org.apache.spark.rdd.MetricsContext._
      sortedRdd.instrument
    }
    instrumentedRdd.mapPartitions(genotypeIterator)
  }

  protected[genotyping] def genotypeSite(region: ReferenceRegion,
                                         referenceObservation: Observation,
                                         alleleObservations: Iterable[AlleleObservation]): Option[VariantContext]

  private def genotypeIterator(iter: Iterator[(ReferencePosition, Observation)]): Iterator[VariantContext] = {
    // do we have any elements in this iterator?
    if (iter.hasNext) {
      // peek at the head element
      val (headPos, headObs) = iter.next

      // and initialize state
      var windowStart = headPos.pos
      var windowEnd = windowStart + headObs.length
      var obsList = List(headObs)

      def flush(): Option[VariantContext] = {
        val region = ReferenceRegion(headPos.referenceName, windowStart, windowEnd)

        // update observations
        val (refObs, obsIterable) = updateObservations(region, obsList)

        // compute genotypes
        genotypeSite(region,
          refObs,
          obsIterable)
      }

      // start processing observation windows!
      iter.flatMap(kv => {
        val (pos, obs) = kv

        // does this observation start past the end of this window? if so, flush our observations
        // and rebuild the window
        if (pos.pos >= windowEnd) {
          // flush the window
          val vc = flush()

          // reinitialize state
          windowStart = pos.pos
          windowEnd = windowStart + obs.length
          obsList = List(obs)

          // emit variant context
          vc
        } else {
          // make possible update to window end
          windowEnd = max(windowEnd, pos.pos + obs.length)

          // append observation to list
          obsList = obs :: obsList

          // nothing to emit
          None
        }
      }) ++ flush().toIterator
    } else {
      Iterator()
    }
  }

  private def updateObservations(region: ReferenceRegion,
                                 obs: List[Observation]): (Observation, Iterable[AlleleObservation]) = {
    // extract references
    val ref = obs.flatMap(o => o match {
      case ao: AlleleObservation => None
      case oo: Observation       => Some(oo)
    }).sortBy(_.pos)
      .map(_.allele)
      .mkString
    val refLen = ref.length
    val refPos = ReferencePosition(region.referenceName, region.start)

    // put read observations together
    val readObservations = obs.flatMap(o => o match {
      case ao: AlleleObservation => Some(ao)
      case oo: Observation       => None
    })

    // if the reference length is greater than 1, we need to expand the values
    val finalObservations = if (refLen > 1) {
      readObservations.groupBy(_.readId)
        .map(ro => {
          // sort the read observations
          val sortedRo = ro._2.sortBy(_.pos)

          // get the start and end position
          val start = sortedRo.head.pos.pos
          val end = sortedRo.last.pos.pos + sortedRo.last.length

          // reduce down to get the sequence
          val sequence = (ref.take((start - region.start).toInt) +
            sortedRo.map(_.allele).mkString +
            ref.takeRight((region.end - end).toInt))
          val phred = sortedRo.map(_.phred).sum / sortedRo.length

          // make observation
          AlleleObservation(refPos,
            refLen,
            sequence,
            phred,
            sortedRo.head.mapq,
            sortedRo.head.onNegativeStrand,
            sortedRo.head.sample,
            sortedRo.head.readId)
        })
    } else {
      readObservations
    }

    // emit reference observation
    (new Observation(refPos, ref), finalObservations)
  }
}
