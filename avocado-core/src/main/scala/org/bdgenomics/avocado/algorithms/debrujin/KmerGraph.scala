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
package org.bdgenomics.avocado.algorithms.debrujin

import org.apache.spark.Logging
import org.bdgenomics.adam.models.{ ReferencePosition, ReferenceRegion }
import org.bdgenomics.avocado.algorithms.hmm.{
  HMMAligner,
  TransitionMatrixConfiguration
}
import org.bdgenomics.avocado.models.{ AlleleObservation, Observation }
import org.bdgenomics.formats.avro.AlignmentRecord
import scala.annotation.tailrec
import scala.collection.mutable.HashMap
import scala.math.abs

object KmerGraph extends Logging {

  /**
   * Creates a new de Brujin Graph that shows connections between k-mers. k-mers are inserted
   * from a read dataset that is provided.
   *
   * @param kmerLength Kmer length.
   * @param readLength Read length.
   * @param regionLength Length of the active region.
   * @param reference Reference sequence.
   * @param reads Reads to insert into graph.
   * @param flankLength Length of flanking sequence to take from reference sequence.
   * @param maxEntries Maximum number of times a loop can be entered.
   * @param lowCoverageTrimmingThreshold Threshold for trimming nodes that are not covered by many reads.
   * @return Returns a new de Brujin graph.
   */
  def apply(kmerLength: Int,
            references: Seq[(ReferenceRegion, String)],
            reads: Seq[AlignmentRecord]): Iterable[KmerGraph] = {

    assert(references.length == 1, "For now, only support a single reference in the graph.")

    log.info("Building indexed de Bruijn graph for region " + references.map(_._1).mkString +
      " from " + reads.size + " reads.")

    @tailrec def addReferenceKmers(iter: Iterator[String],
                                   pos: ReferencePosition,
                                   lastKmer: Kmer,
                                   kmerMap: HashMap[String, Kmer]) {
      // do we have k-mers left?
      if (iter.hasNext) {
        var newKmer: Kmer = null
        val ks = iter.next

        // if we have a predecessor, populate the predecessor fields
        if (lastKmer != null) {
          val kl = List(lastKmer)
          newKmer = Kmer(ks, Some(pos), predecessors = kl)
          lastKmer.successors = newKmer :: lastKmer.successors
        } else {
          newKmer = Kmer(ks, Some(pos))
        }

        // add the kmer to the graph
        kmerMap(ks) = newKmer

        // update the position for the next k-mer
        val newPos = ReferencePosition(pos.referenceName, pos.pos + 1L)

        addReferenceKmers(iter, newPos, newKmer, kmerMap)
      }
    }

    @tailrec def addReadKmers(kmerIter: Iterator[String],
                              qualIter: Iterator[Char],
                              mapq: Int,
                              readId: Long,
                              isNegativeStrand: Boolean,
                              lastKmer: Kmer,
                              kmerMap: HashMap[String, Kmer]) {
      // do we have k-mers left?
      if (kmerIter.hasNext) {
        assert(qualIter.hasNext, "Should still have qualities as long as we've got k-mers.")
        var newKmer: Kmer = null
        val ks = kmerIter.next
        val q = qualIter.next

        // is this a valid k-mer?
        if (!ks.contains('N')) {
          // have we already seen this k-mer?
          if (kmerMap.contains(ks)) {
            newKmer = kmerMap(ks)

            // add phred score and mapq
            newKmer.phred = q.toInt :: newKmer.phred
            newKmer.mapq = mapq :: newKmer.mapq
            newKmer.readId = readId :: newKmer.readId
            newKmer.isNegativeStrand = isNegativeStrand :: newKmer.isNegativeStrand
            assert(newKmer.multiplicity > 0, newKmer.toDetailedString)

            // do we have a predecessor? if so, perform book keeping...
            if (lastKmer != null) {
              // if we have a predecessor, and it isn't in the k-mer map, then add it
              if (newKmer.predecessors.filter(_.kmerSeq == lastKmer.kmerSeq).length == 0) {
                newKmer.predecessors = lastKmer :: newKmer.predecessors
              }

              // if this k-mer isn't in the successor list, then add it
              if (lastKmer.successors.filter(_.kmerSeq == ks).length == 0) {
                lastKmer.successors = newKmer :: lastKmer.successors
              }
            }
          } else {
            val phred = List(q.toInt)
            val mapQ = List(mapq)
            val rid = List(readId)
            val ns = List(isNegativeStrand)

            // if we have a predecessor, populate the predecessor fields
            if (lastKmer != null) {
              val kl = List(lastKmer)
              newKmer = Kmer(ks, phred = phred, mapq = mapQ, readId = rid, isNegativeStrand = ns, predecessors = kl)
              assert(newKmer.multiplicity > 0, newKmer.toDetailedString)
              lastKmer.successors = newKmer :: lastKmer.successors
            } else {
              newKmer = Kmer(ks, phred = phred, mapq = mapQ, readId = rid, isNegativeStrand = ns)
            }

            // add the kmer to the graph
            kmerMap(ks) = newKmer
          }
        }

        addReadKmers(kmerIter, qualIter, mapq, readId, isNegativeStrand, newKmer, kmerMap)
      }
    }

    // group reads by sample
    val readsBySample = reads.groupBy(_.getRecordGroupSample)
    assert(readsBySample.size >= 1, "Must have at least one sample to create a graph.")
    var rId = 0L

    // per sample, construct a debrujin graph
    readsBySample.map(kv => {
      val (sample, sampleReads) = kv

      val kmerMap = HashMap[String, Kmer]()

      // per reference we are passed, add a reference k-mer
      references.foreach(p => {
        val (region, reference) = p

        addReferenceKmers(reference.sliding(kmerLength),
          ReferencePosition(region.referenceName, region.start),
          null,
          kmerMap)
      })

      // loop over reads and collect statistics
      sampleReads.foreach(r => {
        addReadKmers(r.getSequence.toString.sliding(kmerLength),
          r.getQual.toString.toIterator,
          r.getMapq,
          rId,
          r.getReadNegativeStrand,
          null,
          kmerMap)
      })
      rId += 1

      // build the graph for this sample
      new KmerGraph(kmerMap.values.toArray, kmerLength, sample, references)
    })
  }
}

/**
 * Graph showing connections between kmers.
 *
 * @param kmers k-mers which were observed.
 */
class KmerGraph(protected val kmers: Array[Kmer],
                protected val kmerLength: Int,
                val sample: String,
                protected val references: Seq[(ReferenceRegion, String)]) extends Serializable with Logging {

  // get reference start/end positions
  private val starts = references.map(kv => ReferencePosition(kv._1.referenceName, kv._1.start)).toSet
  private val ends = references.map(kv => ReferencePosition(kv._1.referenceName, kv._1.end - kmerLength)).toSet

  // source/sink kmers
  private val allSourceKmers = kmers.filter(_.predecessors.length == 0)
  private val allSinkKmers = kmers.filter(_.successors.length == 0)
  private val sourceKmers = kmers.filter(_.refPos.fold(false)(starts(_)))
  private val sinkKmers = kmers.filter(_.refPos.fold(false)(ends(_)))

  // hmm based aligner for complex variants
  // use default alignment parameters, except infinite padding penalty
  private val hmm = new HMMAligner(new TransitionMatrixConfiguration(LOG_PADDING_PENALTY = Double.NegativeInfinity))

  override def toString(): String = {
    "Sources: " + sourceKmers.map(_.kmerSeq).reduce(_ + ", " + _) + "\n" +
      "Sinks: " + sinkKmers.map(_.kmerSeq).reduce(_ + ", " + _) + "\n" +
      kmers.map(_.toDetailedString).reduce(_ + "\n" + _)
  }

  /**
   * Prints this de Brujin graph out in Graphviz format. This can be used with the Dot software
   * to visualize graph creation.
   *
   * @return Returns a string describing this de Brujin graph as a directed Graphviz graph.
   */
  def toDot(): String = {
    "digraph kg { \n" +
      kmers.map(_.toDot).reduce(_ + "\n" + _) + "\n}\n"
  }

  /**
   * Converts this graph of k-mers into a set of observations. The observations
   * include both the reference threading as well as the actual read derived observations
   * that can be used for genotyping. This is performed without realignment.
   *
   * @return Returns a seq of observations.
   */
  def toObservations: Seq[Observation] = {
    def buildReadObservations(kmer: Kmer,
                              pos: ReferencePosition,
                              length: Int,
                              allele: String,
                              activeReads: Set[Long]): Seq[Observation] = {
      val setSize = if (activeReads == null) {
        0
      } else {
        activeReads.size
      }
      (0 until kmer.multiplicity).flatMap(i => {
        if (activeReads == null || (i < setSize && activeReads(kmer.readId(i)))) {
          Some(AlleleObservation(pos,
            length,
            allele,
            kmer.phred(i),
            kmer.mapq(i),
            kmer.isNegativeStrand(i),
            sample,
            kmer.readId(i)))
        } else {
          None
        }
      })
    }

    def buildReferenceObservations(kmer: Kmer): Seq[Observation] = {
      val site = Seq(new Observation(kmer.refPos.get, kmer.kmerSeq.take(1)))
      // did we have any reads covering this site?
      // if not, just emit that we observed the site
      if (kmer.multiplicity == 0) {
        site
      } else {
        site ++ buildReadObservations(kmer, kmer.refPos.get, 1, kmer.kmerSeq.take(1), null)
      }
    }

    @tailrec def crawl(context: BranchContext,
                       observations: Seq[Observation],
                       branches: List[BranchContext],
                       referenceSites: Set[ReferencePosition]): Seq[Observation] = {

      if (context == null) {
        observations
      } else {
        val (newContext,
          newObservations,
          newBranches,
          newSites) = context match {
          case a: Allele => {
            val newAllele = a.allele + a.kmer.kmerSeq.take(1)
            val newPending = a.kmer :: a.pending

            // is this allele a spur? if so, move on, else continue building the allele
            if (a.kmer.successors.isEmpty) {
              // do we have branches remaining?
              if (branches.isEmpty) {
                (null, observations, branches, referenceSites)
              } else {
                (branches.head, observations, branches.drop(1), referenceSites)
              }
            } else {
              // build next step
              val ctxs = a.kmer.successors.map(nextKmer => {
                if (nextKmer.isReference) {
                  ClosedAllele(nextKmer, newAllele, a.branchPoint, newPending, a.activeReads)
                } else {
                  Allele(nextKmer, newAllele, a.branchPoint, newPending, a.activeReads)
                }
              })

              // return next step
              (ctxs.head,
                observations,
                branches ::: ctxs.drop(1),
                referenceSites)
            }
          }
          case ca: ClosedAllele => {
            // where have we connected back to, and where did we start?
            val refSink = ca.kmer.refPos.get

            // reverse pending k-mers
            val reversed = ca.pending.reverse

            // what is the length of the allele?
            // a snp has length 0, while a n-base indel has length n
            val offset = reversed.length - (refSink.pos - ca.branchPoint.pos).toInt + 1
            val alleleLength = abs(offset)
            lazy val gapLength = if (ca.branchPoint.referenceName == refSink.referenceName) {
              (refSink.pos - ca.branchPoint.pos).toInt - 1
            } else {
              Int.MaxValue
            }

            // do we have an indel/snp or a deletion? deletions must be handeled specially.
            // s/mnp: allele length = 0
            // simple insert: # kmers = gap length + allele length
            // simple deletion: # kmers < kmer length
            if (offset < 0 && reversed.length < kmerLength) {
              // create observations of the alt allele
              val altObs = reversed.takeRight(1)
                .map(buildReadObservations(_,
                  ReferencePosition(refSink.referenceName,
                    refSink.pos - 1),
                  alleleLength,
                  ("_" * alleleLength),
                  ca.activeReads))
                .reduce(_ ++ _)

              // now, count (k - 1) from the start and build reference observations
              var refPoint = ReferencePosition(ca.branchPoint.referenceName,
                ca.branchPoint.pos + 1)
              val refObs = reversed.dropRight(1)
                .map(k => {
                  // take observations
                  val obs = buildReadObservations(k,
                    refPoint,
                    1,
                    k.kmerSeq.take(1),
                    ca.activeReads)

                  // increment position
                  refPoint = ReferencePosition(refPoint.referenceName,
                    refPoint.pos + 1)

                  // emit observations
                  obs
                }).reduce(_ ++ _)

              // combine observations into a new observation set and recurse
              val newObs = observations ++ altObs ++ refObs
              (Reference(ca.kmer),
                newObs,
                branches,
                referenceSites)
            } else if (offset == 0 ||
              (offset > 0 && gapLength + 1 == kmerLength)) {
              // create observations of the alt allele
              val altObs = if (offset == 0) {
                // if we have a s/mnp, we must index into the allele
                var idx = 0
                val alt = ca.allele.drop(kmerLength - alleleLength - 1)
                reversed.drop(kmerLength - 1)
                  .map(kmer => {
                    val obs = buildReadObservations(kmer,
                      ReferencePosition(ca.branchPoint.referenceName,
                        ca.branchPoint.pos + kmerLength - offset + idx),
                      alleleLength,
                      alt(idx).toString,
                      ca.activeReads)

                    // increment index into allele
                    idx += 1

                    obs
                  }).reduce(_ ++ _)
              } else {
                reversed.drop(kmerLength - 1)
                  .map(buildReadObservations(_,
                    ReferencePosition(ca.branchPoint.referenceName,
                      ca.branchPoint.pos + kmerLength - offset),
                    alleleLength,
                    ca.allele.drop(kmerLength - alleleLength - 1),
                    ca.activeReads))
                  .reduce(_ ++ _)
              }

              // now, count (k - 1) from the start and build reference observations
              var refPoint = ReferencePosition(ca.branchPoint.referenceName,
                ca.branchPoint.pos + 1)
              val refObs = reversed.take(kmerLength - 1)
                .map(k => {
                  // take observations
                  val obs = buildReadObservations(k,
                    refPoint,
                    1,
                    k.kmerSeq.take(1),
                    ca.activeReads)

                  // increment position
                  refPoint = ReferencePosition(refPoint.referenceName,
                    refPoint.pos + 1)

                  // emit observations
                  obs
                }).reduce(_ ++ _)

              // combine observations into a new observation set and recurse
              val newObs = observations ++ altObs ++ refObs
              (Reference(ca.kmer),
                newObs,
                branches,
                referenceSites)
            } else {
              // complex alleles represent up in this else statement!
              // processing complex alleles is a 3 step process. we must:
              // 1. recover the reference
              assert(references.length == 1)
              val (refRegion, ref) = references.head
              val refSeq = ref.drop((ca.branchPoint.pos - refRegion.start).toInt + 1).take(gapLength)

              // 2. realign the _allele_ to the reference at that site
              val alignment = hmm.alignSequences(refSeq, ca.allele, null)

              // 3. step through and process the local alignment information

              def kmersToObservations(kmers: Iterator[Kmer],
                                      length: Int,
                                      pos: ReferencePosition): (Seq[Observation], Kmer) = {
                assert(length > 0 && kmers.hasNext)

                // duplicate iterator
                val (k1, k2) = kmers.duplicate
                var kmer: Kmer = null // we must return the last kmer, in case we are followed by a deletion

                // generate allele
                val allele = k1.take(length).map(k => {
                  kmer = k
                  k.kmerSeq.take(1)
                }).reduce(_ + _)

                // generate observations
                (k2.take(length).flatMap(k => buildReadObservations(k,
                  pos,
                  1,
                  allele,
                  ca.activeReads)).toSeq, kmer)
              }

              @tailrec def processComplexAllele(lastKmer: Kmer,
                                                alignmentIter: Iterator[Char],
                                                kmerIter: Iterator[Kmer],
                                                pos: ReferencePosition,
                                                alleleLength: Int = 0,
                                                obs: Seq[Observation] = Seq()): Seq[Observation] = {
                // append the current observations
                // and return the current head kmer, in case we have a deletion
                def addObservations(): (Seq[Observation], Kmer) = {
                  if (alleleLength > 0) {
                    val (newObs, k) = kmersToObservations(kmerIter, alleleLength, pos)
                    (obs ++ newObs, k)
                  } else {
                    (obs, lastKmer)
                  }
                }

                // do we have any more state in the alignment? if not, process the last allele.
                if (!alignmentIter.hasNext) {
                  // when we run out of alignment states, we should either have:
                  // - 1 n-length (n > 0) allele left, and thus n k-mers left
                  // - or, have just processed a 0 length allele (a deletion)
                  assert(kmerIter.hasNext || alleleLength == 0)
                  val finalObs = addObservations()

                  // now, we must be out of k-mers, or something has gone wrong
                  assert(!kmerIter.hasNext)
                  finalObs._1
                } else {
                  // get the next alignment state
                  val alignmentState = alignmentIter.next

                  // if we are here, we must have k-mers left, unless we're in a long deletion
                  assert(kmerIter.hasNext || (alleleLength == 0 && alignmentState == 'D'))

                  // handle the different alignment states:
                  // 1. if we have an alignment match (=/X), we advance the position by 1
                  // 2. if we have a deletion, we emit a length 0 allele and advance the position by 1
                  // 3. if we have an insert, we increment the allele length and don't move
                  val (newLength, newPos, newObs, newKmer) = alignmentState match {
                    case '=' | 'X' => {
                      val (ob, km) = addObservations()
                      (1, ReferencePosition(pos.referenceName, pos.pos + 1), ob, km)
                    }
                    case 'P' | 'D' => {
                      val currPos = ReferencePosition(pos.referenceName, pos.pos + 1)
                      val (ob, km) = addObservations()
                      val delObs = buildReadObservations(km, currPos, 1, "_", ca.activeReads)
                      (0, currPos, ob ++ delObs, km)
                    }
                    case 'I' => {
                      (alleleLength + 1, pos, obs, lastKmer)
                    }
                    case _ => {
                      throw new IllegalArgumentException("Received illegal alignment state in " + alignment)
                    }
                  }

                  // recurse until we are done processing
                  processComplexAllele(newKmer,
                    alignmentIter,
                    kmerIter,
                    newPos,
                    newLength,
                    newObs)
                }
              }

              val newObs = processComplexAllele(reversed.head,
                alignment.alignmentStateSequence.toIterator,
                reversed.toIterator,
                ReferencePosition(ca.branchPoint.referenceName,
                  ca.branchPoint.pos)) ++ observations

              (Reference(ca.kmer),
                newObs,
                branches,
                referenceSites)
            }
          }
          case r: Reference => {
            val pos = r.kmer.refPos.get

            // have we seen this position before?
            // if so, don't process it again
            if (referenceSites(pos)) {
              // do we have remaining sites?
              val (next, newBranches) = if (branches.isEmpty) {
                (null, branches)
              } else {
                (branches.head, branches.drop(1))
              }

              (next,
                observations,
                newBranches,
                referenceSites)
            } else {
              // build observations
              val newObservations = buildReferenceObservations(r.kmer)

              // add new position to set
              val newSet: Set[ReferencePosition] = (referenceSites + pos)

              // if we have a successor, we'll take that and continue on
              // else, we'll move to the next branch
              val (next, newBranches) = if (r.kmer.successors.filter(_.refPos.isDefined).isEmpty && branches.isEmpty) {
                (null, branches)
              } else if (r.kmer.successors.filter(_.refPos.isDefined).isEmpty) {
                (branches.head, branches.drop(1))
              } else {
                (Reference(r.kmer.successors.filter(_.refPos.isDefined).head),
                  r.kmer.successors.filter(_.refPos.isEmpty).map(v => Allele(v, pos)) ::: branches)
              }

              (next,
                observations ++ newObservations,
                newBranches,
                newSet)
            }
          }
        }

        // recurse and compute next iteration
        crawl(newContext,
          newObservations,
          newBranches,
          newSites)
      }
    }

    if (sourceKmers.size > 0) {
      log.info("Started crawling " + references.map(_._1).mkString(", "))
      val obs = crawl(Reference(sourceKmers.head),
        Seq(),
        sourceKmers.drop(1).map(p => Reference(p)).toList,
        Set())
      log.info("Finished crawling " + references.map(_._1).mkString(", "))
      obs
    } else {
      log.warn("No sources seen on " + references.map(_._1).mkString(", "))
      Seq()
    }
  }

  def size: Int = kmers.length

  def nonRefSize: Int = kmers.filter(_.refPos.isEmpty).length

  def sources: Int = sourceKmers.length

  def sinks: Int = sinkKmers.length

  def spurs: Int = allSourceKmers.length + allSinkKmers.length - sources - sinks
}
