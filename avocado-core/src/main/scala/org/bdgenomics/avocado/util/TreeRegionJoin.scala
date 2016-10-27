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
package org.bdgenomics.avocado.util

import com.esotericsoftware.kryo.io.{ Input, Output }
import com.esotericsoftware.kryo.{ Kryo, Serializer }
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{
  ReferenceRegion,
  ReferenceRegionSerializer
}
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.utils.intervaltree.IntervalTree
import scala.annotation.tailrec
import scala.reflect.ClassTag

/**
 * Companion object for building a forest from an RDD.
 */
private[avocado] object Forest extends Serializable {

  /**
   * Sorts the RDD and collects it to build the tree.
   *
   * @param rdd RDD to build a forest from.
   * @return The forest built from this RDD.
   */
  def apply[T: ClassTag](rdd: RDD[(ReferenceRegion, T)]): Forest[T] = BuildingTrees.time {
    val sortedArray = SortingRightSide.time {
      rdd.sortByKey()
        .collect
    }

    Forest(sortedArray)
  }
}

/**
 * Originally, a forest was a collection of trees.
 * Alas, we have no trees anymore.
 * I blame global warming.
 *
 * @param array An array of values for the left side of the join. We require
 *   this array to be sorted.
 */
private[avocado] case class Forest[T: ClassTag](array: Array[(ReferenceRegion, T)]) {

  val length = array.length
  val midpoint = pow2ceil()

  @tailrec private def pow2ceil(i: Int = 1): Int = {
    if (2 * i >= length) {
      i
    } else {
      pow2ceil(2 * i)
    }
  }

  @tailrec private def binarySearch(rr: ReferenceRegion,
                                    idx: Int = 0,
                                    step: Int = midpoint): Option[Int] = {
    if (rr.overlaps(array(idx)._1)) {
      Some(idx)
    } else if (step == 0) {
      None
    } else {
      val stepIdx = idx + step
      val nextIdx = if (stepIdx >= length ||
        (!rr.overlaps(array(stepIdx)._1) &&
          rr.compareTo(array(stepIdx)._1) < 0)) {
        idx
      } else {
        stepIdx
      }
      binarySearch(rr, nextIdx, step / 2)
    }
  }

  @tailrec private def expand(rr: ReferenceRegion,
                              idx: Int,
                              step: Int,
                              list: List[T] = List.empty): List[T] = {
    if (idx < 0 ||
      idx >= length ||
      !rr.overlaps(array(idx)._1)) {
      list
    } else {
      expand(rr, idx + step, step, array(idx)._2 :: list)
    }
  }

  /**
   * @param rr The reference region to grab.
   * @return All keys that overlap the reference region.
   */
  def get(rr: ReferenceRegion): Iterable[T] = {

    val optIdx = binarySearch(rr)

    optIdx.toIterable
      .flatMap(idx => {
        expand(rr, idx, -1) ::: expand(rr, idx + 1, 1)
      })
  }
}

class ForestSerializer[T: ClassTag, TS <: Serializer[T]](
    private val tSerializer: TS) extends Serializer[Forest[T]] {

  private val rrSerializer = new ReferenceRegionSerializer()

  def tTag: ClassTag[T] = implicitly[ClassTag[T]]

  def write(kryo: Kryo, output: Output, obj: Forest[T]) {

    // we will use the array length to allocate an array on read
    output.writeInt(obj.length)

    // loop and write elements
    (0 until obj.length).foreach(idx => {
      rrSerializer.write(kryo, output, obj.array(idx)._1)
      tSerializer.write(kryo, output, obj.array(idx)._2)
    })
  }

  def read(kryo: Kryo, input: Input, klazz: Class[Forest[T]]): Forest[T] = {

    // read the array size and allocate
    val length = input.readInt()
    val array = new Array[(ReferenceRegion, T)](length)

    // loop and read
    (0 until length).foreach(idx => {
      array(idx) = (rrSerializer.read(kryo, input, classOf[ReferenceRegion]),
        tSerializer.read(kryo, input, tTag.runtimeClass.asInstanceOf[Class[T]]))
    })

    Forest[T](array)
  }
}

/**
 * Implements a shuffle free broadcast region join.
 *
 * The broadcast values are stored in a sorted array. It was going to be an
 * ensemble of interval trees, but, that didn't work out.
 */
object TreeRegionJoin extends Serializable {

  /**
   * Performs an inner region join between two RDDs, and groups by the
   * value on the right side of the join.
   *
   * @param leftRdd RDD on the left side of the join. Will be collected to the
   *   driver and broadcast.
   * @param rightRdd RDD on the right side of the join.
   * @return Returns an RDD where each element is a value from the right RDD,
   *   along with all values from the left RDD that it overlapped.
   */
  def joinAndGroupByRight[T: ClassTag, U](
    leftRdd: RDD[(ReferenceRegion, T)],
    rightRdd: RDD[(ReferenceRegion, U)]): RDD[(Iterable[T], U)] = TreeJoin.time {

    // build the tree from the left RDD
    val tree = Forest(leftRdd)

    RunningMapSideJoin.time {
      // broadcast this tree
      val broadcastTree = leftRdd.context
        .broadcast(tree)

      // map and join
      rightRdd.flatMap(kv => {
        val (rr, u) = kv

        // what values keys does this overlap in the tree?
        val overlappingValues = broadcastTree.value
          .get(rr)

        // did we get any overlapping values?
        if (overlappingValues.nonEmpty) {
          Some((overlappingValues, u))
        } else {
          None
        }
      })
    }
  }
}
