package org.bdgenomics.avocado.genotyping

import java.io

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
 * Counts Bases (A C G T) as well as ins / del at each position covered by at least one read in the input alignment file
 */
object PosAlleleStats extends Serializable with Logging {

  case class AlleleCounts(var A: Int, var C: Int, var G: Int, var T: Int, var N: Int, var ins: Int, var del: Int)

  def runPosAlleleStats(alignmentRecordRDD: AlignmentRecordRDD, phredThreshold: Int = 33): RDD[(org.bdgenomics.adam.models.ReferencePosition, org.bdgenomics.avocado.genotyping.PosAlleleStats.AlleleCounts)] = {

    def posAlleleStatsMapPartitions(in: Iterator[AlignmentRecord]): Iterator[(ReferencePosition, AlleleCounts)] = {

      // hash of allele counters per genomics position
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
          val sequence: String = read.getSequence
          val qual: String = read.getQual

          println(sequence)
          println(qual)

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

          def incrementBaseCount(i: Char, pos: ReferencePosition, qual: Int): Unit = {

            if (!x.isDefinedAt(pos)) {
              x(pos) = AlleleCounts(0, 0, 0, 0, 0, 0, 0)
            }

            if (qual >= phredThreshold) {
              if (i.toUpper == 'A') {
                x(pos).A += 1
              } else if (i.toUpper == 'C') {
                x(pos).C += 1
              } else if (i.toUpper == 'G') {
                x(pos).G += 1
              } else if (i.toUpper == 'T') {
                x(pos).T += 1
              }
            } else {
              x(pos).N += 1
            }
          }

          @tailrec def recordAlleles(iter: Iterator[ObservationOperator],
                                     lastRef: String = ""): Unit = {

            if (!iter.hasNext) {
              None
            } else {
              // pop from the iterator
              val obs = iter.next
              // update the list and position and advance
              val (nextRef) = obs match {
                case Clipped(_, _) => {
                  (lastRef)
                }
                case Match(length, optRef) => {
                  val kv = optRef.fold({
                    (0 until length).foreach(i => {
                      incrementBaseCount(sequence(idx + i), ReferencePosition(contigName, pos + i), qual(idx + i).toInt - 33)
                      //println("test qual i:" + i + " idx:" + idx + " alt: " + qual(i).toString + " fixed_alt: " + qual(i + idx) + " qual: " + (qual(idx + i).toInt - 33))
                    })
                    sequence(idx + length - 1).toString
                  })((ref: String) => {
                    (0 until length).foreach(i => {
                      //println("test qual i:" + i + " idx:" + idx + " ref:" + ref(i).toString() + " alt: " + qual(i).toString + " fixed_alt: " + qual(i + idx) + " qual: " + (qual(idx + i).toInt - 33))
                      incrementBaseCount(sequence(idx + i), ReferencePosition(contigName, pos + i), qual(idx + i).toInt - 33)
                    })
                    ref.last.toString
                  })
                  pos += length
                  idx += length
                  kv
                }
                case Insertion(length) => {
                  val insQuals = qual.substring(idx - 1, idx + length).map(_.toInt - 33).sum / length
                  if (insQuals >= phredThreshold) { x(ReferencePosition(contigName, pos - 1)).ins += 1 }
                  idx += length
                  (lastRef)
                }
                case Deletion(ref) => {
                  val delLength = ref.size
                  if (qual(idx - 1).toInt - 33 >= phredThreshold) { x(ReferencePosition(contigName, pos)).del += 1 }
                  pos += delLength
                  (ref.last.toString)
                }
              }

              // recurse
              recordAlleles(iter, nextRef)
            }
          }

          recordAlleles(opsIter)
        }

      })

      x.toIterator
    }

    alignmentRecordRDD.rdd.mapPartitions(posAlleleStatsMapPartitions)
      .reduceByKey((x, y) => AlleleCounts(x.A + y.A,
        x.C + y.C,
        x.G + y.G,
        x.T + y.T,
        x.N + y.N,
        x.ins + y.ins,
        x.del + y.del))
  }

}

