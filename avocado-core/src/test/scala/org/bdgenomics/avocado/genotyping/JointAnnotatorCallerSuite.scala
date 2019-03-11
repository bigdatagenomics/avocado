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

import org.bdgenomics.adam.models.VariantContext
import org.bdgenomics.avocado.AvocadoFunSuite
import org.bdgenomics.formats.avro.{
  Genotype,
  GenotypeAllele,
  Variant,
  VariantCallingAnnotations
}
import org.bdgenomics.utils.misc.MathUtils
import scala.collection.JavaConversions._

class JointAnnotatorCallerSuite extends AvocadoFunSuite {

  val baseGt = Genotype.newBuilder
    .setReferenceName("chr1")
    .setStart(1000)
    .setEnd(1001)
    .setVariant(Variant.newBuilder
      .setAlternateAllele("A")
      .setReferenceAllele("G")
      .build)
    .build

  test("discard reference site") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.REF))
        .build)

    val optNewSite = JointAnnotatorCaller.annotateSite(
      VariantContext.buildFromGenotypes(genotypes))

    assert(optNewSite.isEmpty)
  }

  test("calculate MAF for all called genotypes") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .build,
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.REF))
        .build)
    val af = JointAnnotatorCaller.calculateMinorAlleleFrequency(
      VariantContext.buildFromGenotypes(genotypes))

    assert(MathUtils.fpEquals(af, 0.25))
  }

  test("calculate MAF ignoring uncalled genotypes") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.REF,
          GenotypeAllele.REF,
          GenotypeAllele.ALT))
        .build,
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.OTHER_ALT))
        .build,
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.NO_CALL, GenotypeAllele.NO_CALL))
        .build)
    val af = JointAnnotatorCaller.calculateMinorAlleleFrequency(
      VariantContext.buildFromGenotypes(genotypes))

    assert(MathUtils.fpEquals(af, 0.25))
  }

  test("roll up variant annotations from a single genotype") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setReadDepth(10)
        .setReferenceReadDepth(5)
        .setAlternateReadDepth(5)
        .setStrandBiasComponents(Seq(4, 1, 2, 3)
          .map(i => i: java.lang.Integer))
        .build)

    val annotation = JointAnnotatorCaller.calculateAnnotations(
      VariantContext.buildFromGenotypes(genotypes))

    assert(annotation.getReadDepth === 10)
    assert(annotation.getReferenceReadDepth === 5)
    assert(annotation.getReferenceForwardReadDepth === 4)
    assert(annotation.getReferenceReverseReadDepth === 1)
    assert(annotation.getForwardReadDepth === 6)
    assert(annotation.getReverseReadDepth === 4)
  }

  test("roll up variant annotations across multiple genotypes") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setReadDepth(10)
        .setReferenceReadDepth(5)
        .setAlternateReadDepth(5)
        .setStrandBiasComponents(Seq(4, 1, 2, 3)
          .map(i => i: java.lang.Integer))
        .build,
      Genotype.newBuilder(baseGt)
        .setReadDepth(7)
        .setReferenceReadDepth(5)
        .setAlternateReadDepth(2)
        .setStrandBiasComponents(Seq(3, 2, 1, 1)
          .map(i => i: java.lang.Integer))
        .build,
      Genotype.newBuilder(baseGt)
        .setReadDepth(9)
        .setReferenceReadDepth(9)
        .setAlternateReadDepth(0)
        .build)

    val annotation = JointAnnotatorCaller.calculateAnnotations(
      VariantContext.buildFromGenotypes(genotypes))

    assert(annotation.getReadDepth === 26)
    assert(annotation.getReferenceReadDepth === 19)
    assert(annotation.getReferenceForwardReadDepth === 7)
    assert(annotation.getReferenceReverseReadDepth === 3)
    assert(annotation.getForwardReadDepth === 10)
    assert(annotation.getReverseReadDepth === 7)
  }

  test("recalling genotypes is a no-op for no calls and complex hets") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.NO_CALL))
        .build,
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.OTHER_ALT, GenotypeAllele.NO_CALL))
        .build)

    val newGenotypes = JointAnnotatorCaller.recallGenotypes(
      VariantContext.buildFromGenotypes(genotypes), 0.5)

    genotypes.zip(newGenotypes)
      .foreach(p => assert(p._1 === p._2))
  }

  test("recall a genotype so that the state changes") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .setGenotypeLikelihoods(Seq(0.0, 0.0, 0.0)
          .map(d => d: java.lang.Double))
        .build)

    val newGenotypes = JointAnnotatorCaller.recallGenotypes(
      VariantContext.buildFromGenotypes(genotypes), 0.9)

    assert(newGenotypes.size === 1)
    val newGenotype = newGenotypes.head
    assert(newGenotype.getAlleles.count(_ == GenotypeAllele.ALT) === 2)
    val posteriors = newGenotype.getVariantCallingAnnotations
      .getGenotypePosteriors
    assert(posteriors.size === 3)
    assert(MathUtils.fpEquals(posteriors.get(0).toDouble, -4.60517018599))
    assert(MathUtils.fpEquals(posteriors.get(1).toDouble, -1.71479842809))
    assert(MathUtils.fpEquals(posteriors.get(2).toDouble, -0.21072103131))
    val priors = newGenotype.getVariantCallingAnnotations
      .getGenotypePriors
    assert(priors.size === 3)
    assert(MathUtils.fpEquals(priors.get(0).toDouble, -4.60517018599))
    assert(MathUtils.fpEquals(priors.get(1).toDouble, -1.71479842809))
    assert(MathUtils.fpEquals(priors.get(2).toDouble, -0.21072103131))
    assert(newGenotype.getGenotypeQuality === 6)
  }

  test("allele frequency being outside of (0.0, 1.0) just computes posteriors") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setAlleles(Seq(GenotypeAllele.REF, GenotypeAllele.ALT))
        .setGenotypeLikelihoods(Seq(0.0, 1.0, 0.0)
          .map(d => d: java.lang.Double))
        .build)

    val newGenotypes = JointAnnotatorCaller.recallGenotypes(
      VariantContext.buildFromGenotypes(genotypes), 0.0)

    assert(newGenotypes.size === 1)
    val newGenotype = newGenotypes.head
    assert(newGenotype.getAlleles.count(_ == GenotypeAllele.ALT) === 1)
    val posteriors = newGenotype.getVariantCallingAnnotations
      .getGenotypePosteriors
    assert(posteriors.size === 3)
    assert(MathUtils.fpEquals(posteriors.get(0).toDouble, -1.55144471393))
    assert(MathUtils.fpEquals(posteriors.get(1).toDouble, -0.55144471393))
    assert(MathUtils.fpEquals(posteriors.get(2).toDouble, -1.55144471393))
    assert(newGenotype.getVariantCallingAnnotations
      .getGenotypePriors
      .isEmpty)
  }

  test("compute variant quality from a single genotype") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setVariantCallingAnnotations(
          VariantCallingAnnotations.newBuilder
            .setGenotypePosteriors(Seq(-4.60517018599,
              -2.40794560865,
              -0.10536051565)
              .map(_.toFloat)
              .map(f => f: java.lang.Float))
            .build)
        .build)

    assert(MathUtils.fpEquals(JointAnnotatorCaller.computeQuality(genotypes),
      20.0))
  }

  test("compute variant quality from multiple genotypes") {
    val genotypes = Seq(
      Genotype.newBuilder(baseGt)
        .setVariantCallingAnnotations(
          VariantCallingAnnotations.newBuilder
            .setGenotypePosteriors(Seq(-4.60517018599,
              -2.40794560865,
              -0.10536051565)
              .map(_.toFloat)
              .map(f => f: java.lang.Float))
            .build)
        .build,
      Genotype.newBuilder(baseGt)
        .setVariantCallingAnnotations(
          VariantCallingAnnotations.newBuilder
            .setGenotypePosteriors(Seq(-6.90775527898,
              -0.00200200267
                - 6.90775527898)
              .map(_.toFloat)
              .map(f => f: java.lang.Float))
            .build)
        .build)

    assert(MathUtils.fpEquals(JointAnnotatorCaller.computeQuality(genotypes),
      50.0))
  }
}
