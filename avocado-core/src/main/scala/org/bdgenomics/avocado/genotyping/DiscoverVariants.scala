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
import org.apache.spark.sql.SQLContext
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variant.VariantRDD
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.models.{
  Clipped,
  Deletion,
  Insertion,
  Match,
  ObservationOperator
}
import org.bdgenomics.formats.avro.{ AlignmentRecord, Variant }
import org.bdgenomics.utils.misc.Logging
import scala.annotation.tailrec

/**
 * Discovers the variants present in a set of aligned reads.
 *
 * Useful for force-calling variants.
 */
object DiscoverVariants extends Serializable with Logging {

  /**
   * Discovers all variants in an RDD of reads.
   *
   * @param rdd RDD of reads.
   * @param optPhredThreshold An optional threshold that discards all variants
   *   not supported by bases of at least a given phred score.
   * @return Returns an RDD of variants.
   */
  private[avocado] def apply(
    aRdd: AlignmentRecordRDD,
    optPhredThreshold: Option[Int] = None,
    optMinObservations: Option[Int] = None): VariantRDD = DiscoveringVariants.time {

    VariantRDD(variantsInRdd(aRdd.rdd,
      optPhredThreshold = optPhredThreshold,
      optMinObservations = optMinObservations),
      aRdd.sequences)
  }

  /**
   * Discovers all variants in an RDD of reads.
   *
   * @param rdd RDD of reads.
   * @param optPhredThreshold An optional threshold that discards all variants
   *   not supported by bases of at least a given phred score.
   * @return Returns an RDD of variants.
   */
  private[genotyping] def variantsInRdd(
    rdd: RDD[AlignmentRecord],
    optPhredThreshold: Option[Int] = None,
    optMinObservations: Option[Int] = None): RDD[Variant] = {

    // if phred threshold is unset, set to 0
    val phredThreshold = optPhredThreshold.getOrElse(0)

    val variantRdd = rdd.flatMap(variantsInRead(_, phredThreshold))

    // convert to dataframe
    val sqlContext = SQLContext.getOrCreate(rdd.context)
    import sqlContext.implicits._
    val variantDs = sqlContext.createDataFrame(variantRdd)

    // count by variant and remove
    val uniqueVariants = optMinObservations.fold({
      variantDs.distinct
    })(mo => {
      variantDs.groupBy(variantDs("contigName"),
        variantDs("start"),
        variantDs("end"),
        variantDs("referenceAllele"),
        variantDs("alternateAllele"))
        .count()
        .where($"count" > mo)
        .drop("count")
    })

    uniqueVariants.as[DiscoveredVariant]
      .rdd
      .map(_.toVariant)
  }

  /**
   * Discovers the variants in a single read.
   *
   * @param read Aligned read to look for variants in.
   * @param phredThreshold A threshold that discards all variants not supported
   *   by bases of at least a given phred score.
   * @return Returns a collection containing all the variants in a read.
   */
  private[genotyping] def variantsInRead(read: AlignmentRecord,
                                         phredThreshold: Int): Iterable[DiscoveredVariant] = {

    if (!read.getReadMapped) {
      Iterable.empty
    } else {
      // extract the alignment blocks from the read
      val ops = try {
        ObservationOperator.extractAlignmentOperators(read)
      } catch {
        case t: Throwable => {
          log.warn("Extracting alignment operators from %s failed with %s.".format(
            read.getReadName, t))
          Iterable.empty
        }
      }

      // where are we on the reference and in the read?
      var pos = read.getStart.toInt
      var idx = 0

      // get the read sequence, contig, etc
      val sequence = read.getSequence
      val qual = read.getQual
      val contigName = read.getContigName

      // advance to the first alignment match
      @tailrec def fastForward(
        iter: BufferedIterator[ObservationOperator]): Iterator[ObservationOperator] = {

        if (!iter.hasNext) {
          Iterator()
        } else {
          val stop = iter.head match {
            case Clipped(_, false) => false
            case Clipped(length, true) => {
              idx += length
              false
            }
            case Insertion(length) => {
              idx += length
              false
            }
            case Deletion(ref) => {
              pos += ref.length
              false
            }
            case Match(_, _) => {
              true
            }
          }

          if (stop) {
            iter.toIterator
          } else {
            // pop from the iterator and recurse
            iter.next
            fastForward(iter)
          }
        }
      }

      val opsIter = fastForward(ops.toIterator.buffered)

      // emit variants
      @tailrec def emitVariants(
        iter: Iterator[ObservationOperator],
        lastRef: String = "",
        variants: List[DiscoveredVariant] = List.empty): Iterable[DiscoveredVariant] = {

        if (!iter.hasNext) {
          variants.toIterable
        } else {

          // pop from the iterator
          val obs = iter.next

          // update the list and position and advance
          val (nextRef, nextVariants) = obs match {
            case Clipped(_, _) => {
              (lastRef, variants)
            }
            case Match(length, optRef) => {
              val kv = optRef.fold({
                (sequence(idx + length - 1).toString, variants)
              })(ref => {
                val newVars = (0 until length).flatMap(i => {
                  if (qual(i).toInt - 33 >= phredThreshold) {
                    Some(DiscoveredVariant(
                      contigName,
                      pos + i,
                      pos + i + 1,
                      ref(i).toString,
                      sequence(idx + i).toString))
                  } else {
                    None
                  }
                }).toList ::: variants
                (ref.last.toString, newVars)
              })
              pos += length
              idx += length
              kv
            }
            case Insertion(length) => {
              val insQuals = qual.substring(idx - 1, idx + length).map(_.toInt - 33).sum / length
              val newVar = if (insQuals >= phredThreshold) {
                DiscoveredVariant(
                  contigName,
                  pos - 1,
                  pos,
                  lastRef,
                  sequence.substring(idx - 1, idx + length)) :: variants
              } else {
                variants
              }
              idx += length
              (lastRef, newVar)
            }
            case Deletion(ref) => {
              val delLength = ref.size
              val newVar = if (qual(idx - 1).toInt - 33 >= phredThreshold) {
                DiscoveredVariant(
                  contigName,
                  pos - 1,
                  pos + delLength,
                  lastRef + ref,
                  sequence.substring(idx - 1, idx)) :: variants
              } else {
                variants
              }
              pos += delLength
              (ref.last.toString, newVar)
            }
          }

          // recurse
          emitVariants(iter, nextRef, nextVariants)
        }
      }

      emitVariants(opsIter)
    }
  }
}
