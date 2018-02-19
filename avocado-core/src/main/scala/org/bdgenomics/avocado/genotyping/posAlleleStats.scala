package org.bdgenomics.avocado.genotyping

import org.apache.spark.rdd.RDD
import org.apache.spark.sql.SQLContext
import org.bdgenomics.adam.models.ReferencePosition
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.adam.rdd.variant.VariantRDD
import org.bdgenomics.avocado.Timers._
import org.bdgenomics.avocado.models.{ Clipped, Deletion, Insertion, Match, ObservationOperator }
import org.bdgenomics.formats.avro.{ AlignmentRecord, Variant }
import org.bdgenomics.utils.misc.Logging

import scala.annotation.tailrec
import scala.collection.mutable

/**
 * Discovers the variants present in a set of aligned reads.
 *
 * Useful for force-calling variants.
 */
object PosAlleleStats extends Serializable with Logging {

  case class AlleleCounts(var A: Int, var C: Int, var G: Int, var T: Int)

  def posAlleleStatsMapPartitions(in: Iterator[AlignmentRecord]): Iterator[(ReferencePosition, AlleleCounts)] = {

    //  val x: mutable.Map[ReferencePosition, mutable.HashMap[Char,Int]] = new scala.collection.mutable.HashMap[org.bdgenomics.adam.models.ReferencePosition, mutable.HashMap[Char,Int]]
    val x: mutable.Map[ReferencePosition, AlleleCounts] = new scala.collection.mutable.HashMap[org.bdgenomics.adam.models.ReferencePosition, AlleleCounts]

    in.foreach(read => {
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

        def incrementBaseCount(i: Char, pos: ReferencePosition): Unit = {

          if (!x.isDefinedAt(pos)) {
            x(pos) = AlleleCounts(0, 0, 0, 0)
          }

          if (i.toUpper == 'A') { x(pos).A += 1 }
          else if (i.toUpper == 'C') { x(pos).C += 1 }
          else if (i.toUpper == 'C') { x(pos).G += 1 }
          else if (i.toUpper == 'T') { x(pos).G += 1 }
        }

        @tailrec def recordAlleles(iter: Iterator[ObservationOperator],
                                   lastRef: String = "",
                                   variants: List[DiscoveredVariant] = List.empty,
                                   phredThreshold: Int = 30): Unit = {

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
                  (0 until length).foreach(i => if (qual(i).toInt - 33 >= phredThreshold) {

                    incrementBaseCount(sequence(idx + i), ReferencePosition(contigName, pos + i))
                  })

                  //x(ReferencePosition(contigName,pos+1))=(1,1,1,1) })
                  /*
                  (0 until length).flatMap( i => if (qual(i).toInt - 33 >= phredThreshold) {
                    x(ReferencePosition(contigName,pos+1))=(1,1,1,1)
                    None } else {
                    None} )
                  */

                  (sequence(idx + length - 1).toString, variants)
                })((ref: String) => {
                  val newVars = (0 until length).flatMap(i => {
                    if (qual(i).toInt - 33 >= phredThreshold) {
                      ///x(ReferencePosition(contigName,pos+i)) = (1,1,1,1)
                      incrementBaseCount(sequence(idx + i), ReferencePosition(contigName, pos + i))
                      Some(DiscoveredVariant(
                        contigName,
                        pos + i,
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
            recordAlleles(iter, nextRef, nextVariants)
          }
        }

        recordAlleles(opsIter)
      }

    }) //foreach read

    //val z: Iterator[(ReferencePosition, AlleleCounts)] = x.toIterator
    x.toIterator
  }

}