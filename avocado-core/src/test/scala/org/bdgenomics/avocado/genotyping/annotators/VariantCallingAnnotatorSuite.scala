package org.bdgenomics.avocado.genotyping.annotators

import org.bdgenomics.avocado.algorithms.math.MathTestUtils
import org.bdgenomics.avocado.genotyping.BiallelicGenotyper
import org.scalatest.FunSuite
import org.bdgenomics.formats.avro.{ Variant, VariantCallingAnnotations }
import org.bdgenomics.adam.models.{ SequenceRecord, SequenceDictionary, ReferencePosition }
import org.bdgenomics.avocado.algorithms.math.LogToPhred.log2phred
import org.bdgenomics.avocado.models.AlleleObservation

import scala.math.{ log, pow, sqrt }

class VariantCallingAnnotatorSuite extends FunSuite {
  val ba = new BiallelicGenotyper(SequenceDictionary(SequenceRecord("ctg", 1000L)),
    Map(("ctg", 1000L)))
  val ref = "A"
  val alt = "G"
  val variant = Variant.newBuilder()
    .setReferenceAllele(ref)
    .setAlternateAllele(alt)
    .build()

  val variantQuality = 35d

  def createAllele(base: String, mapQ: Int, strand: Boolean): AlleleObservation = {
    AlleleObservation(ReferencePosition("ctg", 0L),
      1,
      base,
      30,
      Some(mapQ),
      strand,
      true,
      1,
      "mySample",
      1L)
  }

  test("annotate variant quality by depth") {
    val observed = Iterable(createAllele(ref, 20, true),
      createAllele(ref, 30, true),
      createAllele(alt, 40, true))

    val likelihoods = ba.scoreGenotypeLikelihoods(ref, alt, observed)._2
    val annotator = VariantCallingAnnotator(variant, observed, variantQuality, likelihoods)
    MathTestUtils.assertAlmostEqual(annotator.variantQualityByDepth(), variantQuality)
  }

  test("annotate read depth") {
    val observed = Iterable(createAllele(ref, 20, true),
      createAllele(ref, 30, true),
      createAllele(alt, 40, true))
    val likelihoods = ba.scoreGenotypeLikelihoods(ref, alt, observed)._2
    val annotator = VariantCallingAnnotator(variant, observed, variantQuality, likelihoods)
    assert(annotator.readDepth() == 3)
  }

  test("annotate map quality rank sum") {
    val refBaseQ = List(20, 25, 26, 30, 32, 40, 47, 50, 53, 60)
    val altBaseQ = List(0, 7, 10, 17, 20, 21, 30, 34, 40, 45)
    val observed = refBaseQ.map(createAllele(ref, _, true)) ++ altBaseQ.map(createAllele(alt, _, false))
    val likelihoods = ba.scoreGenotypeLikelihoods(ref, alt, observed)._2
    val annotator = VariantCallingAnnotator(variant, observed, variantQuality, likelihoods)
    val rs = annotator.mapQRankSum()
    assert(rs.isDefined)
    MathTestUtils.assertAlmostEqual(rs.get, -2.154f, 1e-3)
  }

  test("can't calculate map quality rank sum for empty ref") {
    val refBaseQ = List()
    val altBaseQ = List(0, 7, 10, 17, 20, 21, 30, 34, 40, 45)
    val observed = refBaseQ.map(createAllele(ref, _, true)) ++ altBaseQ.map(createAllele(alt, _, false))
    val likelihoods = ba.scoreGenotypeLikelihoods(ref, alt, observed)._2
    val annotator = VariantCallingAnnotator(variant, observed, variantQuality, likelihoods)
    val rs = annotator.mapQRankSum()
    assert(rs.isEmpty)
  }

  test("can't calculate map quality rank sum for empty alt") {
    val refBaseQ = List(20, 25, 26, 30, 32, 40, 47, 50, 53, 60)
    val altBaseQ = List()
    val observed = refBaseQ.map(createAllele(ref, _, true)) ++ altBaseQ.map(createAllele(alt, _, false))
    val likelihoods = ba.scoreGenotypeLikelihoods(ref, alt, observed)._2
    val annotator = VariantCallingAnnotator(variant, observed, variantQuality, likelihoods)
    val rs = annotator.mapQRankSum()
    assert(rs.isEmpty)
  }

  test("annotate fisher strand bias value") {
    val observed = Iterable(createAllele(ref, 20, true),
      createAllele(ref, 30, true),
      createAllele(alt, 40, false),
      createAllele(alt, 40, false))
    val likelihoods = ba.scoreGenotypeLikelihoods(ref, alt, observed)._2
    val annotator = VariantCallingAnnotator(variant, observed, variantQuality, likelihoods)
    MathTestUtils.assertAlmostEqual(annotator.fisherStrandBiasValue(), log2phred(log(1 / 6f)))
  }

  test("annotate rms map q") {
    val observed = Iterable(createAllele(ref, 20, true),
      createAllele(ref, 30, true),
      createAllele(alt, 40, false),
      createAllele(alt, 40, false))
    val likelihoods = ba.scoreGenotypeLikelihoods(ref, alt, observed)._2
    val annotator = VariantCallingAnnotator(variant, observed, variantQuality, likelihoods)
    val rmsMapQ = sqrt((pow(20d, 2d) + pow(30d, 2d) + pow(40d, 2d) + pow(40d, 2d)) / 4d)
    MathTestUtils.assertAlmostEqual(annotator.rmsMapQ(), rmsMapQ)
  }
}
