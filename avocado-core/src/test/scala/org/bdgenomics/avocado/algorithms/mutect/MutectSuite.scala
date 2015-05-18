package org.bdgenomics.avocado.algorithms.mutect

import org.bdgenomics.adam.models.ReferencePosition
import org.bdgenomics.avocado.models.AlleleObservation
import org.scalatest.FunSuite

/**
 * Created by john on 5/17/15.
 */
class MutectSuite extends FunSuite {
  val pos = ReferencePosition("ctg", 101l)
  val mt = new Mutect()
  //Demo data
  val all_c = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = "C",
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = true,
      sample = "tumor",
      readId = 1L)
  }

  val some_muts = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = true,
      sample = "tumor",
      readId = 1L)
  }

  test("Sites with no mutations should not have any variants returned") {
    val result = mt.detect(pos, "C", all_c)
    assert(result.isEmpty, "Result should be empty")
  }

  test("Sites with mutations should return a variant") {
    val result = mt.detect(pos, "C", some_muts)
    assert(result.isDefined, "Result should not be empty")
    val (variant, genotype) = result.get
    assert(variant.getAlternateAllele === "A", "Alternate allele should be a")
    assert(variant.getReferenceAllele === "C", "Reference allele should be c")
    assert(variant.getContig.getContigName === "ctg", "Contig should be 'ctg'")
    assert(variant.getStart === 101l, "Variant should be at position 101")
    assert(variant.getEnd === 102l, "Variant should be at position 101")
  }

}
