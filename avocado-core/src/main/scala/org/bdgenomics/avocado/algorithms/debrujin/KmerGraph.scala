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
import org.bdgenomics.avocado.Timers._
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

    var observationEstimate = 0

    @tailrec def addReferenceKmers(iter: Iterator[String],
                                   pos: ReferencePosition,
                                   lastKmer: Kmer,
                                   kmerMap: HashMap[String, Kmer]) {
      // do we have k-mers left?
      if (iter.hasNext) {
        var newKmer: Kmer = null
        val ks = iter.next

        // increment our observation pointer
        observationEstimate += 1

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

        // increment our observation pointer
        observationEstimate += 1

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

      // reset observation count
      observationEstimate = 0

      val kmerMap = HashMap[String, Kmer]()

      BuildingReferenceGraph.time {
        // per reference we are passed, add a reference k-mer
        references.foreach(p => {
          val (region, reference) = p

          addReferenceKmers(reference.sliding(kmerLength),
            ReferencePosition(region.referenceName, region.start),
            null,
            kmerMap)
        })
      }

      AddingReadsToGraph.time {
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
      }

      ConstructingGraph.time {
        // build the graph for this sample
        new KmerGraph(kmerMap.values.toArray,
          kmerLength,
          sample,
          references,
          Some(observationEstimate))
      }
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
                protected val references: Seq[(ReferenceRegion, String)],
                hintedCoverage: Option[Int] = None) extends Serializable with Logging {

  // get reference start/end positions
  private val starts = references.map(kv => ReferencePosition(kv._1.referenceName, kv._1.start)).toSet
  private val ends = references.map(kv => ReferencePosition(kv._1.referenceName, kv._1.end - kmerLength)).toSet

  // source/sink kmers
  private val allSourceKmers = kmers.filter(_.predecessors.length == 0)
  private val allSinkKmers = kmers.filter(_.successors.length == 0)
  private val sourceKmers = kmers.filter(_.refPos.fold(false)(starts(_)))
  private val sinkKmers = kmers.filter(_.refPos.fold(false)(ends(_)))

  // how many observations will we have in this region?
  private lazy val coverage = hintedCoverage.getOrElse(kmers.map(k => {
    k.refPos.fold(0)(v => 1) + k.multiplicity
  }).sum)

  // a summary string for this region
  private lazy val refString = references.map(_._1).mkString(", ")

  override def toString(): String = {
    "Sources: " + sourceKmers.map(_.kmerSeq).mkString(", ") + "\n" +
      "Sinks: " + sinkKmers.map(_.kmerSeq).mkString(", ") + "\n" +
      kmers.map(_.toDetailedString).mkString("\n")
  }

  /**
   * Prints this de Brujin graph out in Graphviz format. This can be used with the Dot software
   * to visualize graph creation.
   *
   * @return Returns a string describing this de Brujin graph as a directed Graphviz graph.
   */
  def toDot(): String = {
    "digraph kg { \n" +
      kmers.map(_.toDot).mkString("\n") + "\n}\n"
  }

  /**
   * Converts this graph of k-mers into a set of observations. The observations
   * include both the reference threading as well as the actual read derived observations
   * that can be used for genotyping. This is performed without realignment.
   *
   * @return Returns a seq of observations.
   */
  def toObservations: Seq[Observation] = {
    // preallocate an array for observations
    var observations = new Array[Observation](coverage)
    var obsIdx = 0
    var obsArrayLen = coverage

    def updateObservations(os: Iterable[Observation]): Unit = UpdatingObservations.time {
      if (obsIdx + os.size >= obsArrayLen) {
        ExtendingArray.time {
          val increment = obsArrayLen >> 2
          log.warn("Extending observation array in " + refString + " by " + increment +
            " from " + obsArrayLen + " to " + (obsArrayLen += increment))
          observations = observations.padTo(obsArrayLen, null)
        }
      }

      // loop and add observations
      os.foreach(o => {
        observations(obsIdx) = o
        obsIdx += 1
      })
    }

    def buildReadObservations(kmer: Kmer,
                              pos: ReferencePosition,
                              length: Int,
                              allele: String,
                              activeReads: Set[Long]): Seq[AlleleObservation] = ReadObservations.time {
      val setSize = if (activeReads == null) {
        0
      } else {
        activeReads.size
      }
      var seen = 0
      (0 until kmer.multiplicity).flatMap(i => {
        if (activeReads == null || (seen < setSize && activeReads(kmer.readId(i)))) {
          seen += 1
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

    def buildReferenceObservations(kmer: Kmer): Seq[Observation] = ReferenceObservations.time {
      val site = Seq(new Observation(kmer.refPos.get, kmer.kmerSeq.take(1)))
      // did we have any reads covering this site?
      // if not, just emit that we observed the site
      if (kmer.multiplicity == 0) {
        site
      } else {
        site ++ (buildReadObservations(kmer, kmer.refPos.get, 1, kmer.kmerSeq.take(1), null)
          .map(_.asInstanceOf[Observation]))
      }
    }

    @tailrec def crawl(context: BranchContext,
                       branches: List[BranchContext],
                       referenceSites: Set[ReferencePosition]) {

      if (context != null) {
        val (newContext,
          newBranches,
          newSites) = Stepping.time {
          context match {
            case a: Allele => CrawlAllele.time {
              val newAllele = a.allele + a.kmer.kmerSeq.take(1)
              val newPending = a.kmer :: a.pending

              // is this allele a spur? if so, move on, else continue building the allele
              if (a.kmer.successors.isEmpty) {
                // do we have branches remaining?
                if (branches.isEmpty) {
                  (null, branches, referenceSites)
                } else {
                  (branches.head, branches.drop(1), referenceSites)
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
                  branches ::: ctxs.drop(1),
                  referenceSites)
              }
            }
            case ca: ClosedAllele => ClosingAllele.time {
              // where have we connected back to, and where did we start?
              val refSink = ca.kmer.refPos.get

              // reverse pending k-mers
              val reversed = ca.pending.reverse

              // what is the length of the allele?
              // a snp has length 0, while a n-base indel has length n
              val refGap = (refSink.pos - ca.branchPoint.pos).toInt
              val offset = reversed.length - refGap + 1
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
                val delBase = reversed.head.kmerSeq.dropRight(1).takeRight(1)

                // create observations of the alt allele
                val altObs = reversed.takeRight(1)
                  .flatMap(buildReadObservations(_,
                    ReferencePosition(refSink.referenceName,
                      refSink.pos - 1),
                    alleleLength + 1,
                    (delBase),
                    ca.activeReads))
                  .map(_.asInstanceOf[Observation])

                // now, count (k - 1) from the start and build reference observations
                var refPoint = ReferencePosition(ca.branchPoint.referenceName,
                  ca.branchPoint.pos + 1)
                val refObs = reversed.dropRight(1)
                  .flatMap(k => {
                    // take observations
                    val obs = buildReadObservations(k,
                      refPoint,
                      1,
                      k.kmerSeq.take(1),
                      ca.activeReads)
                      .map(_.asInstanceOf[Observation])

                    // increment position
                    refPoint = ReferencePosition(refPoint.referenceName,
                      refPoint.pos + 1)

                    // emit observations
                    obs
                  })

                // combine observations into a new observation set and recurse
                val newObs = altObs ++ refObs
                updateObservations(newObs)
                (Reference(ca.kmer),
                  branches,
                  referenceSites)
              } else {
                // create observations of the alt allele
                val altObs = if (offset == 0) {
                  // if we have a s/mnp, we must index into the allele
                  var idx = 0
                  val alt = ca.allele.drop(kmerLength - alleleLength - 1)
                  reversed.drop(kmerLength - 1)
                    .flatMap(kmer => {
                      val obs = buildReadObservations(kmer,
                        ReferencePosition(ca.branchPoint.referenceName,
                          ca.branchPoint.pos + kmerLength - offset + idx),
                        1,
                        alt(idx).toString,
                        ca.activeReads)
                        .map(_.asInstanceOf[Observation])

                      // increment index into allele
                      idx += 1

                      obs
                    })
                } else {
                  // are we in an insert, or a complex allele? if we are in an insert, we have
                  // allele length of 1, if we're in a complex allele, we emit the short
                  // haplotype corresponding to the complex allele
                  val (altLength, altAllele, rp) = if (offset > 0 && gapLength + 1 == kmerLength) {
                    // a simple insert left aligns to the reference position prior to the insert
                    (1,
                      ca.allele.drop(kmerLength - alleleLength - 1),
                      ReferencePosition(ca.branchPoint.referenceName,
                        ca.branchPoint.pos + kmerLength - offset))
                  } else {
                    // a complex allele covers multiple reference positions
                    (abs(offset - 1),
                      ca.allele.drop(kmerLength - 1),
                      ReferencePosition(ca.branchPoint.referenceName,
                        ca.branchPoint.pos + kmerLength))
                  }
                  reversed.drop(kmerLength - 1)
                    .flatMap(buildReadObservations(_,
                      rp,
                      altLength,
                      altAllele,
                      ca.activeReads))
                    .groupBy(ao => ao.readId)
                    .map(kv => {
                      val (_, observations) = kv

                      // seed from the head observation
                      val headObs = observations.head

                      // we need to average out the phred scores
                      val numObservations = observations.length
                      val phred = observations.map(_.phred).sum

                      // create new observation that is the average of the others
                      AlleleObservation(headObs.pos,
                        headObs.length,
                        headObs.allele,
                        phred / numObservations,
                        headObs.mapq,
                        headObs.onNegativeStrand,
                        headObs.sample,
                        headObs.readId).asInstanceOf[Observation]
                    })
                }

                // now, count (k - 1) from the start and build reference observations
                var refPoint = ReferencePosition(ca.branchPoint.referenceName,
                  ca.branchPoint.pos + 1)
                val refObs = reversed.take(kmerLength - 1)
                  .flatMap(k => {
                    // take observations
                    val obs = buildReadObservations(k,
                      refPoint,
                      1,
                      k.kmerSeq.take(1),
                      ca.activeReads)
                      .map(_.asInstanceOf[Observation])

                    // increment position
                    refPoint = ReferencePosition(refPoint.referenceName,
                      refPoint.pos + 1)

                    // emit observations
                    obs
                  })

                // combine observations into a new observation set and recurse
                val newObs = altObs ++ refObs
                updateObservations(newObs)
                (Reference(ca.kmer),
                  branches,
                  referenceSites)
              }
            }
            case r: Reference => CrawlReference.time {
              val pos = r.kmer.refPos.get

              // have we seen this position before?
              // if so, don't process it again
              if (CheckingIfHaveSeen.time { referenceSites(pos) }) {
                SiteSeenBefore.time {
                  // do we have remaining sites?
                  val (next, newBranches) = if (branches.isEmpty) {
                    (null, branches)
                  } else {
                    (branches.head, branches.drop(1))
                  }

                  (next,
                    newBranches,
                    referenceSites)
                }
              } else {
                ProcessingUnseenSite.time {
                  // build observations
                  val newObservations = ProcessingObservations.time {
                    buildReferenceObservations(r.kmer)
                  }

                  // add new position to set
                  val newSet: Set[ReferencePosition] = RebuildingSet.time { (referenceSites + pos) }

                  // if we have a successor, we'll take that and continue on
                  // else, we'll move to the next branch
                  val (next, newBranches) = PickingBranch.time {
                    if (r.kmer.successors.filter(_.refPos.isDefined).isEmpty && branches.isEmpty) {
                      (null, branches)
                    } else if (r.kmer.successors.filter(_.refPos.isDefined).isEmpty) {
                      (branches.head, branches.drop(1))
                    } else {
                      (Reference(r.kmer.successors.filter(_.refPos.isDefined).head),
                        r.kmer.successors.filter(_.refPos.isEmpty).map(v => Allele(v, pos)) ::: branches)
                    }
                  }

                  // update observations and return
                  updateObservations(newObservations)
                  BuildingBranchInfo.time {
                    (next,
                      newBranches,
                      EvaluatingSet.time { newSet })
                  }
                }
              }
            }
          }
        }

        // recurse and compute next iteration
        crawl(newContext,
          newBranches,
          newSites)
      }
    }

    if (sourceKmers.size > 0) {
      log.info("Started crawling " + refString)
      crawl(Reference(sourceKmers.head),
        sourceKmers.drop(1).map(p => Reference(p)).toList,
        Set())
      log.info("Finished crawling " + refString)
      observations.take(obsIdx).toSeq
    } else {
      log.warn("No sources seen on " + refString)
      Seq()
    }
  }

  def size: Int = kmers.length

  def nonRefSize: Int = kmers.filter(_.refPos.isEmpty).length

  def sources: Int = sourceKmers.length

  def sinks: Int = sinkKmers.length

  def spurs: Int = allSourceKmers.length + allSinkKmers.length - sources - sinks
}
