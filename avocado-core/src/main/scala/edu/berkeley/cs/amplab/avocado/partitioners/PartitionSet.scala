/*
 * Copyright (c) 2014. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.berkeley.cs.amplab.avocado.partitioners

import edu.berkeley.cs.amplab.adam.models.ReferenceRegion
import scala.annotation.tailrec
import scala.collection.immutable.{SortedMap, TreeSet, SortedSet}

class PartitionSet (protected val regionMapping: SortedMap[ReferenceRegion, Int]) extends Serializable {

  /**
   * Merges a region into a sorted set of regions.
   *
   * @note Assumes that the region being merged in is ordered after all regions in the
   * current set of regions.
   * 
   * @param set Sorted set of previously merged regions.
   * @param region Region to merge in to set.
   * @return New set with new region merged in.
   */
  private def mergeRegions (set: TreeSet[ReferenceRegion],
                            region: ReferenceRegion): TreeSet[ReferenceRegion] = {
    if (set.isEmpty) {
      set + region
    } else {
      val t = set.last
      
      if (t.overlaps(region) || t.isAdjacent(region)) {
        set.dropRight(1) + t.merge(region)
      } else {
        set + region
      }
    }
  }

  /**
   * Merges overlapping/adjacent regions together across the set of regions.
   *
   * @param regions Input set of regions.
   * @return Sorted set with overlapping regions merged together.
   */
  private def merge (regions: SortedSet[ReferenceRegion]): TreeSet[ReferenceRegion] = {
    regions.foldLeft(TreeSet[ReferenceRegion]())(mergeRegions)
  }

  lazy val mergedPartitions = merge(regionMapping.keySet)

  /**
   * Returns a list of all integer partition mappings that a region overlaps with.
   *
   * @note Can be overridden if a more performant partition mapping function can be provided.
   * 
   * @param region Region of interest.
   * @return List of all partition indexes that this region overlaps with.
   */
  def getPartition (region: ReferenceRegion): List[Int] = {
    regionMapping.filterKeys(_.overlaps(region))
      .values
      .toList
  }

  /**
   * Returns whether a region is contained in a partition inside of this set.
   *
   * @param region Region of interest.
   * @return True if region overlaps with any region inside of this partition.
   */
  final def isInSet (region: ReferenceRegion): Boolean = {
    !mergedPartitions.filter(_.refId == region.refId)
      .forall(r => !(r.contains(region) || r.overlaps(region)))
  }

  /**
   * Returns whether a region is not wholly contained inside of this set.
   *
   * @param region Region of interest.
   * @return True if region is not wholly contained inside of this set.
   */
  final def isOutsideOfSet (region: ReferenceRegion): Boolean = {
    mergedPartitions.filter(_.refId == region.refId)
      .forall(r => !r.contains(region))
  }

}
