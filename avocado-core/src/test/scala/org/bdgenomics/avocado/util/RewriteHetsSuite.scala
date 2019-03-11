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

import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.adam.models.SequenceDictionary
import org.bdgenomics.adam.rdd.variant.GenotypeDataset
import org.bdgenomics.formats.avro.{ Genotype, GenotypeAllele, Variant }
import scala.collection.JavaConversions._

case class TestRewriteHetsArgs(
    var maxHetSnpAltAllelicFraction: Float = 0.75f,
    var maxHetIndelAltAllelicFraction: Float = 0.65f,
    var disableHetSnpRewriting: Boolean = false,
    var disableHetIndelRewriting: Boolean = false) extends RewriteHetsArgs {
}

class RewriteHetsSuite extends AvocadoFunSuite {

  def buildGt(ref: String,
              alt: String,
              gq: Int,
              dp: Int,
              adp: Int,
              alleles: List[GenotypeAllele]): Genotype = {
    Genotype.newBuilder
      .setGenotypeQuality(gq)
      .setVariant(Variant.newBuilder
        .setReferenceAllele(ref)
        .setAlternateAllele(alt)
        .build)
      .setAlleles(alleles)
      .setReadDepth(dp)
      .setAlternateReadDepth(adp)
      .build
  }

  val goodHetSnp = buildGt("A", "T", 30, 30, 15,
    List(GenotypeAllele.REF, GenotypeAllele.ALT))
  val badHetSnp = buildGt("A", "T", 30, 30, 25,
    List(GenotypeAllele.REF, GenotypeAllele.ALT))
  val goodHetIndel = buildGt("AA", "T", 30, 50, 30,
    List(GenotypeAllele.REF, GenotypeAllele.ALT))
  val badHetIndel = buildGt("A", "TCG", 30, 20, 20,
    List(GenotypeAllele.REF, GenotypeAllele.ALT, GenotypeAllele.ALT))
  val homRefSnp = buildGt("A", "T", 30, 30, 3,
    List(GenotypeAllele.REF, GenotypeAllele.REF))
  val homRefIndel = buildGt("A", "TCC", 30, 50, 5,
    List(GenotypeAllele.REF))
  val homAltSnp = buildGt("A", "T", 30, 30, 29,
    List(GenotypeAllele.ALT))
  val homAltIndel = buildGt("A", "TCC", 30, 50, 45,
    List(GenotypeAllele.ALT, GenotypeAllele.ALT))

  def shouldRewrite(gt: Genotype,
                    snps: Boolean = true,
                    indels: Boolean = true): Boolean = {
    RewriteHets.shouldRewrite(gt, 0.75f, 0.65f, snps, indels)
  }

  test("should rewrite a bad het snp") {
    assert(shouldRewrite(badHetSnp))
  }

  test("should not rewrite het snp if snp filtering is disabled") {
    assert(!shouldRewrite(badHetSnp, snps = false))
  }

  test("should rewrite a bad het indel") {
    assert(shouldRewrite(badHetIndel))
  }

  test("should not rewrite het indel if indel filtering is disabled") {
    assert(!shouldRewrite(badHetIndel, indels = false))
  }

  test("don't rewrite good het calls") {
    assert(!shouldRewrite(goodHetSnp))
    assert(!shouldRewrite(goodHetIndel))
  }

  test("don't rewrite homozygous calls") {
    assert(!shouldRewrite(homRefSnp))
    assert(!shouldRewrite(homRefIndel))
    assert(!shouldRewrite(homAltSnp))
    assert(!shouldRewrite(homAltIndel))
  }

  test("rewrite a het call as a hom alt snp") {
    val rewrittenGt = RewriteHets.rewriteGenotype(badHetSnp)
    assert(rewrittenGt.getVariant.getReferenceAllele === "A")
    assert(rewrittenGt.getVariant.getAlternateAllele === "T")
    assert(Option(rewrittenGt.getGenotypeQuality).isEmpty)
    assert(rewrittenGt.getAlleles.length === 2)
    assert(rewrittenGt.getAlleles.get(0) === GenotypeAllele.ALT)
    assert(rewrittenGt.getAlleles.get(1) === GenotypeAllele.ALT)
    assert(rewrittenGt.getReadDepth === 30)
    assert(rewrittenGt.getAlternateReadDepth === 25)
  }

  def processGenotype(gt: Genotype,
                      snps: Boolean = true,
                      indels: Boolean = true): Genotype = {
    RewriteHets.processGenotype(gt, 0.75f, 0.65f, snps, indels)
  }

  test("processing a valid call should not change the call") {
    val goodCalls = Seq(goodHetSnp, goodHetIndel,
      homRefSnp, homRefIndel,
      homAltSnp, homAltIndel)
    goodCalls.foreach(gt => assert(gt === processGenotype(gt)))
  }

  test("if processing is disabled, don't rewrite bad calls") {
    val badCalls = Seq(badHetSnp, badHetIndel)
    badCalls.foreach(gt => assert(gt === processGenotype(gt,
      snps = false,
      indels = false)))
  }

  test("process a bad het snp call") {
    val rewrittenGt = processGenotype(badHetSnp)
    assert(rewrittenGt != badHetSnp)
    assert(rewrittenGt.getVariant.getReferenceAllele === "A")
    assert(rewrittenGt.getVariant.getAlternateAllele === "T")
    assert(Option(rewrittenGt.getGenotypeQuality).isEmpty)
    assert(rewrittenGt.getAlleles.length === 2)
    assert(rewrittenGt.getAlleles.get(0) === GenotypeAllele.ALT)
    assert(rewrittenGt.getAlleles.get(1) === GenotypeAllele.ALT)
    assert(rewrittenGt.getReadDepth === 30)
    assert(rewrittenGt.getAlternateReadDepth === 25)
  }

  test("process a bad het indel call") {
    val rewrittenGt = processGenotype(badHetIndel)
    assert(rewrittenGt != badHetIndel)
    assert(rewrittenGt.getVariant.getReferenceAllele === "A")
    assert(rewrittenGt.getVariant.getAlternateAllele === "TCG")
    assert(Option(rewrittenGt.getGenotypeQuality).isEmpty)
    assert(rewrittenGt.getAlleles.length === 3)
    assert(rewrittenGt.getAlleles.get(0) === GenotypeAllele.ALT)
    assert(rewrittenGt.getAlleles.get(1) === GenotypeAllele.ALT)
    assert(rewrittenGt.getAlleles.get(2) === GenotypeAllele.ALT)
    assert(rewrittenGt.getReadDepth === 20)
    assert(rewrittenGt.getAlternateReadDepth === 20)
  }

  val genotypes = Seq(badHetSnp, badHetIndel,
    goodHetSnp, goodHetIndel,
    homRefSnp, homRefIndel,
    homAltSnp, homAltIndel)

  def gtRdd: GenotypeDataset = {
    val rdd = sc.parallelize(genotypes)
    GenotypeDataset(rdd,
      SequenceDictionary.empty,
      Seq.empty,
      Seq.empty)
  }

  sparkTest("disable processing for a whole rdd") {
    val rewrittenRdd = RewriteHets(gtRdd,
      TestRewriteHetsArgs(disableHetSnpRewriting = true,
        disableHetIndelRewriting = true))
    val newGts = rewrittenRdd.rdd.collect
    val oldGtSet = genotypes.toSet
    newGts.foreach(gt => assert(oldGtSet(gt)))
  }

  sparkTest("process a whole rdd") {
    val rewrittenRdd = RewriteHets(gtRdd, TestRewriteHetsArgs())
    val newGts = rewrittenRdd.rdd.collect
    val oldGtSet = genotypes.toSet
    val (touchedGts, untouchedGts) = newGts.partition(_.getGenotypeQuality == null)
    assert(untouchedGts.size === 6)
    assert(touchedGts.size === 2)
    untouchedGts.foreach(gt => assert(oldGtSet(gt)))
    touchedGts.foreach(gt => assert(gt.getAlleles.forall(_ == GenotypeAllele.ALT)))
  }
}
