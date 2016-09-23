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
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variation.VariantRDD
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
   * @return Returns an RDD of variants.
   */
  private[avocado] def apply(aRdd: AlignmentRecordRDD): VariantRDD = {

    VariantRDD(variantsInRdd(aRdd.rdd),
      aRdd.sequences)
  }

  /**
   * Discovers all variants in an RDD of reads.
   *
   * @param rdd RDD of reads.
   * @return Returns an RDD of variants.
   */
  private[genotyping] def variantsInRdd(rdd: RDD[AlignmentRecord]): RDD[Variant] = {

    rdd.flatMap(variantsInRead)
      .distinct
  }

  /**
   * Discovers the variants in a single read.
   *
   * @param read Aligned read to look for variants in.
   * @return Returns a collection containing all the variants in a read.
   */
  private[genotyping] def variantsInRead(read: AlignmentRecord): Iterable[Variant] = {

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
      var pos = read.getStart
      var idx = 0

      // get the read sequence, contig, etc
      val sequence = read.getSequence
      val contigName = read.getContigName

      // advance to the first alignment match
      @tailrec def fastForward(
        iter: BufferedIterator[ObservationOperator]): Iterator[ObservationOperator] = {

        if (!iter.hasNext) {
          Iterator()
        } else {
          val stop = iter.head match {
            case Clipped(length, _) => {
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
      @tailrec def emitVariants(iter: Iterator[ObservationOperator],
                                lastRef: String = "",
                                variants: List[Variant] = List.empty): Iterable[Variant] = {

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
                (ref.last.toString,
                  Variant.newBuilder
                  .setContigName(contigName)
                  .setStart(pos)
                  .setEnd(pos + length.toLong)
                  .setReferenceAllele(ref)
                  .setAlternateAllele(sequence.substring(idx, idx + length))
                  .build :: variants)
              })
              pos += length
              idx += length
              kv
            }
            case Insertion(length) => {
              val insVariant = Variant.newBuilder
                .setContigName(contigName)
                .setStart(pos - 1L)
                .setEnd(pos)
                .setReferenceAllele(lastRef)
                .setAlternateAllele(lastRef + sequence.substring(idx, idx + length))
                .build
              idx += length
              (lastRef, insVariant :: variants)
            }
            case Deletion(ref) => {
              val delLength = ref.size
              val delVariant = Variant.newBuilder
                .setContigName(contigName)
                .setStart(pos - 1L)
                .setEnd(pos + delLength.toLong)
                .setReferenceAllele(lastRef + ref)
                .setAlternateAllele(lastRef)
                .build
              pos += delLength
              (ref.last.toString, delVariant :: variants)
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
