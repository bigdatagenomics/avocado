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
package org.bdgenomics.avocado.algorithms.join

import org.apache.spark.{ Logging, Partitioner, RangePartitioner, SparkContext }
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.avocado.algorithms.join.JoinOrderings._
import scala.annotation.tailrec
import scala.collection.mutable.Queue
import scala.reflect.ClassTag
import scala.math.{ max, min }

object ShuffleMultiJoin extends Serializable with Logging {

  private[join] def processPartition[T, U](iterT: Iterator[(ReferenceRegion, T)],
                                           iterU: Iterator[(ReferenceRegion, U)]): Iterator[(T, Iterable[U])] = {
    var haveSeenT = false
    var tQueue: Queue[(ReferenceRegion, T)] = Queue.empty
    var uQueue: Queue[(ReferenceRegion, U)] = Queue.empty

    @tailrec def check(region: ReferenceRegion,
                       rList: List[(T, Iterable[U])] = List.empty): Iterable[(T, Iterable[U])] = {
      if (tQueue.isEmpty || region.overlaps(tQueue.head._1)) {
        rList.toIterable
      } else if (region.referenceName == tQueue.head._1.referenceName &&
        region.start >= tQueue.head._1.end) {
        val (r, t) = tQueue.dequeue
        val newList = (t, uQueue.filter(_._1.overlaps(r))
          .map(_._2)
          .toIterable) :: rList

        // trim from the head of the queue
        if (!tQueue.isEmpty) {
          val newR = tQueue.head._1
          uQueue = uQueue.dropWhile(!_._1.overlaps(newR))
        }

        check(region, newList)
      } else {
        Iterable.empty
      }
    }

    def flush(): Iterator[(T, Iterable[U])] = {
      (tQueue.toIterator ++ iterT)
        .map(kv => {
          val (r, t) = kv
          (t, uQueue.filter(_._1.overlaps(r))
            .map(_._2)
            .toIterable)
        })
    }

    iterU.flatMap(kv => {
      uQueue += kv
      val result = check(kv._1)
      if (tQueue.isEmpty && iterT.hasNext) {
        tQueue += iterT.next
      }
      result
    }) ++ flush()
  }

  def partitionAndMultiJoin[T, U](leftRDD: RDD[(ReferenceRegion, T)],
                                  rightRDD: RDD[(ReferenceRegion, U)])(implicit tManifest: ClassTag[T],
                                                                       uManifest: ClassTag[U]): RDD[(T, Iterable[U])] = {
    // build an index from the right side RDD
    val index = RegionIndexer.buildBySamplingRDD(rightRDD)

    def prepare[V](rdd: RDD[(ReferenceRegion, V)])(implicit vTag: ClassTag[V]): RDD[(ReferenceRegion, V)] = {
      rdd.flatMap(kv => {
        val idxs = index.getPartitions(kv._1)
        if (idxs.isEmpty) {
          log.warn("Dropping %s, which has no valid partition in the join.".format(kv))
        }
        idxs.map(i => ((i, kv._1), kv._2))
      }).repartitionAndSortWithinPartitions(IndexPartitioner(index.numPartitions))
        .map(kv => (kv._1._2, kv._2))
    }

    // flat map over both right and left RDDs
    val preparedLeftRDD = prepare(leftRDD)
    val preparedRightRDD = prepare(rightRDD)

    // then, map partitions
    preparedLeftRDD.zipPartitions(preparedRightRDD)(processPartition)
  }
}
