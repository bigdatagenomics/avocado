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

import org.bdgenomics.adam.models.SequenceDictionary
import org.bdgenomics.adam.rdd.read.AlignmentRecordRDD
import org.bdgenomics.formats.avro.AlignmentRecord

private[avocado] trait PrefilterReadsArgs extends Serializable {

  /**
   * True if a genome build is not from the Genome Reference Consortium.
   */
  var isNotGrc: Boolean

  /**
   * True if we want to restrict our reads to the autosomal chromosomes.
   */
  var autosomalOnly: Boolean

  /**
   * True if we want to keep reads aligned to the mitochondrial chromosome.
   */
  var keepMitochondrialChromosome: Boolean

  /**
   * True if we want to keep reads marked as sequencing duplicates.
   */
  var keepDuplicates: Boolean

  /**
   * The minimum mapping quality of read to keep.
   */
  var minMappingQuality: Int

  /**
   * If true, keeps secondary/supplemental alignments.
   */
  var keepNonPrimary: Boolean
}

/**
 * Reifies an input AlignmentRecordRDD down to the contigs and reads we
 * want to genotype.
 */
private[avocado] object PrefilterReads extends Serializable {

  /**
   * Filters out reads and contigs that should not be processed.
   *
   * @param rdd RDD of reads and associated metadata.
   * @param args Arguments specifying the filters to apply.
   * @return Returns a new AlignmentRecordRDD where reads that we don't want
   *   to use in genotyping have been discarded, and where contigs that we
   *   don't want to genotype have been removed.
   */
  def apply(rdd: AlignmentRecordRDD,
            args: PrefilterReadsArgs): AlignmentRecordRDD = {

    // get filter functions
    val contigFn = contigFilterFn(args)
    val readFn = readFilterFn(args, contigFn)

    // filter contigs and construct a new sequence dictionary
    val sequences = new SequenceDictionary(rdd.sequences
      .records
      .filter(r => contigFn(r.name)))

    // filter reads and construct a new rdd
    rdd.transform(_.filter(readFn))
      .replaceSequences(sequences)
  }

  /**
   * @param args The arguments specifying which contigs to keep.
   * @return Returns a function that returns true if a contig with a given name
   *   should be kept.
   */
  protected[util] def contigFilterFn(args: PrefilterReadsArgs): (String => Boolean) = {
    val fns = if (args.isNotGrc) {
      Iterable(filterNonGrcAutosome(_), filterNonGrcSex(_), filterNonGrcMitochondrial(_))
    } else {
      Iterable(filterGrcAutosome(_), filterGrcSex(_), filterGrcMitochondrial(_))
    }
    val filteredFns = Iterable(true, !args.autosomalOnly, args.keepMitochondrialChromosome)
      .zip(fns)
      .filter(_._1)
      .map(_._2)

    assert(filteredFns.nonEmpty)

    def filterFn(s: String): Boolean = {
      filteredFns.map(fn => fn(s)).reduce(_ || _)
    }

    filterFn(_)
  }

  /**
   * @param args The arguments specifying which reads to keep.
   * @param contigFilterFn A function that determines which contigs should be
   *   kept, given the contig name.
   * @return Returns a function that returns true if a read should be kept.
   */
  protected[util] def readFilterFn(
    args: PrefilterReadsArgs,
    contigFilterFn: (String => Boolean)): (AlignmentRecord => Boolean) = {

    def baseFilterFn(r: AlignmentRecord): Boolean = {
      (filterMapped(r, args.keepNonPrimary) &&
        filterMappingQuality(r, args.minMappingQuality) &&
        contigFilterFn(r.getContigName))
    }

    if (args.keepDuplicates) {
      baseFilterFn(_)
    } else {
      def filterFn(r: AlignmentRecord): Boolean = {
        filterUnique(r) && baseFilterFn(r)
      }

      filterFn(_)
    }
  }

  /**
   * @param read Read to test for filtration.
   * @return Returns true if the read is not a duplicate read.
   */
  protected[util] def filterUnique(read: AlignmentRecord): Boolean = {
    !read.getDuplicateRead
  }

  /**
   * @param read Read to test for filtration.
   * @return Returns true if the read is aligned.
   */
  protected[util] def filterMapped(read: AlignmentRecord,
                                   keepNonPrimary: Boolean): Boolean = {
    read.getReadMapped && (keepNonPrimary || read.getPrimaryAlignment)
  }

  /**
   * @param read Read to test for filtration.
   * @param minMappingQuality Only keep reads with a mapping quality above this
   *   value.
   * @return Returns true if the read is aligned.
   */
  protected[util] def filterMappingQuality(read: AlignmentRecord,
                                           minMappingQuality: Int): Boolean = {
    // if mapq is not set, ignore
    if (read.getMapq == null) {
      true
    } else {
      read.getMapq > minMappingQuality
    }
  }

  /**
   * @param contigName Contig name to test for filtration.
   * @return Returns true if the contig matches the naming scheme for GRCh
   *   autosomal chromosomes.
   */
  protected[util] def filterGrcAutosome(contigName: String): Boolean = {
    contigName != null &&
      contigName.size >= 4 &&
      contigName.startsWith("chr") && contigName.drop(3).forall(_.isDigit)
  }

  /**
   * @param contigName Contig name to test for filtration.
   * @return Returns true if the contig matches the naming scheme for GRCh
   *   sex chromosomes.
   */
  protected[util] def filterGrcSex(contigName: String): Boolean = {
    if (contigName != null &&
      contigName.length == 4 &&
      contigName.startsWith("chr")) {
      contigName(3) == 'X' || contigName(3) == 'Y' ||
        contigName(3) == 'Z' || contigName(3) == 'W'
    } else {
      false
    }
  }

  /**
   * @param contigName Contig name to test for filtration.
   * @return Returns true if the contig matches the GRCh mitochondrial
   *   chromosome name.
   */
  protected[util] def filterGrcMitochondrial(contigName: String): Boolean = {
    contigName != null && contigName == "chrM"
  }

  /**
   * @param contigName Contig name to test for filtration.
   * @return Returns true if the contig matches the naming scheme for HG/UCSC
   *   autosomal chromosomes.
   */
  protected[util] def filterNonGrcAutosome(contigName: String): Boolean = {
    contigName != null && contigName.forall(_.isDigit)
  }

  /**
   * @param contigName Contig name to test for filtration.
   * @return Returns true if the contig matches the naming scheme for HG/UCSC
   *   sex chromosomes.
   */
  protected[util] def filterNonGrcSex(contigName: String): Boolean = {
    contigName != null &&
      (contigName == "X" || contigName == "Y" ||
        contigName == "Z" || contigName == "W")
  }

  /**
   * @param contigName Contig name to test for filtration.
   * @return Returns true if the contig matches the HG/UCSC mitochondrial
   *   chromosome name.
   */
  protected[util] def filterNonGrcMitochondrial(contigName: String): Boolean = {
    contigName != null && contigName == "MT"
  }
}
