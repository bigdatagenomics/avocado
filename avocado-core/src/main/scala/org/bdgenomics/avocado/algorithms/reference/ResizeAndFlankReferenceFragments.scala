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
package org.bdgenomics.avocado.algorithms.reference

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{ ReferenceRegion, SequenceDictionary }
import org.bdgenomics.formats.avro.NucleotideContigFragment

object ResizeAndFlankReferenceFragments extends Serializable {

  def apply(rdd: RDD[NucleotideContigFragment],
            sd: SequenceDictionary,
            fragmentLength: Int,
            flankSize: Int): RDD[NucleotideContigFragment] = {
    rdd.keyBy(ctg => ReferenceRegion(ctg).get)
      .repartitionAndSortWithinPartitions(ReferencePartitioner(sd))
      .mapPartitions(overlapAndResize(_, fragmentLength, flankSize))
  }

  private[reference] def resize(contig: (ReferenceRegion, NucleotideContigFragment),
                                targetLength: Double): Array[(ReferenceRegion, NucleotideContigFragment)] = {
    val (region, fragment) = contig

    // we can't guarantee that we'll split the fragment perfectly into the target length,
    // but we'll pick the closest integer splits
    val numSplits = (region.length.toDouble / targetLength).toInt
    val splitLength = (region.length / numSplits).toInt

    // state for splitter
    val array = new Array[(ReferenceRegion, NucleotideContigFragment)](numSplits)
    var fragmentSequence = fragment.getFragmentSequence
    var start = region.start

    // loop and process splits
    (0 until numSplits).foreach(i => {
      // the last split must take all remaining bases
      val (fragmentRegion, fragmentString) = if (i == numSplits - 1) {
        (ReferenceRegion(region.referenceName, start, start + fragmentSequence.length.toLong),
          fragmentSequence)
      } else {
        val kv = (ReferenceRegion(region.referenceName, start, start + splitLength.toLong),
          fragmentSequence.take(splitLength))

        // update state
        start += splitLength
        fragmentSequence = fragmentSequence.drop(splitLength)

        kv
      }

      // add to the array
      array(i) = (fragmentRegion, NucleotideContigFragment.newBuilder()
        .setFragmentSequence(fragmentString)
        .setFragmentStartPosition(fragmentRegion.start)
        .setContig(fragment.getContig)
        .build())
    })

    // return stuffed array
    array
  }

  private[reference] def overlapAndResize(iter: Iterator[(ReferenceRegion, NucleotideContigFragment)],
                                          fragmentLength: Int,
                                          flankSize: Int): Iterator[NucleotideContigFragment] = {
    // we need to have at least one element in the iterator
    if (iter.hasNext) {
      // first, we split the fragments down into their new size
      val splitFragments = iter.flatMap(resize(_, fragmentLength))

      // now, we apply a window and flank adjacent segments
      var lastFragment = splitFragments.next
      splitFragments.map(f => {
        // grab temp copy; we will overwrite later
        val copyLastFragment = lastFragment

        // are the two fragments adjacent? if so, we must add the flanking sequences
        if (copyLastFragment._1.isAdjacent(f._1)) {
          val lastSequence = copyLastFragment._2.getFragmentSequence
          val currSequence = f._2.getFragmentSequence

          // update fragments with flanking sequences
          copyLastFragment._2.setFragmentSequence(lastSequence + currSequence.take(flankSize))
          copyLastFragment._2.setDescription(Option(copyLastFragment._2.getDescription)
            .fold("rr")(_ + "rr"))
          f._2.setFragmentSequence(lastSequence.takeRight(flankSize) + currSequence)
          f._2.setDescription("f")

          // we must change the start position of the fragment we are appending in front of
          f._2.setFragmentStartPosition(f._2.getFragmentStartPosition - flankSize.toLong)
        }

        // overwrite last fragment
        lastFragment = f

        // emit updated last fragment
        copyLastFragment._2
      }) ++ Iterator(lastFragment._2)
    } else {
      Iterator()
    }
  }
}
