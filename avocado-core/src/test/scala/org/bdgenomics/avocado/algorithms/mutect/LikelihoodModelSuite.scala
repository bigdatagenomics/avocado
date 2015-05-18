package org.bdgenomics.avocado.algorithms.mutect

import org.bdgenomics.adam.models.ReferencePosition
import org.bdgenomics.avocado.algorithms.math.MathTestUtils
import org.bdgenomics.avocado.models.AlleleObservation
import org.scalatest.FunSuite

/**
 * Created by john on 5/17/15.
 */
class LikelihoodModelSuite extends FunSuite {
  //Demo data
  val all_c = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = "C",
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = true,
      sample = "shouldntmatter",
      readId = 1L)
  }

  val some_muts = for (read_id <- 1 to 10) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id <= 3) "A" else "C", //Mutants = A x 3 w/ scores (30,31,32)
      phred = 30 + (read_id - 1), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = true,
      sample = "shouldntmatter",
      readId = 1L)
  }

  val normals_het = for (read_id <- 11 to 20) yield {
    AlleleObservation(pos = ReferencePosition("ctg", 0L),
      length = 1,
      allele = if (read_id % 2 == 0) "A" else "C", //phreads == 31, 33, ... 39
      phred = 30 + (read_id - 11), // 30 to 39
      mapq = Some(30),
      onNegativeStrand = true,
      sample = "shouldntmatter",
      readId = 1L)
  }

  test("Likelihood of simple model, no errors or mutants, all reference sites") {
    val m0likelihood = M0Model.logLikelihood("C", "A", all_c, None)
    val mflikelihood = MfmModel.logLikelihood("C", "A", all_c, None)
    // Assuming f = 0, mflikelihood or m0likelihood
    // in R:
    // options(digits=10)
    // p = seq(30,39)
    // e = 10^(-p/10)
    // log10(prod(1-e))
    // > -0.001901013984
    MathTestUtils.assertAlmostEqual(m0likelihood, -0.001901013984)
    MathTestUtils.assertAlmostEqual(mflikelihood, m0likelihood)

    val mhlikelihood = MHModel.logLikelihood("C", "A", all_c, None)
    // assuming f = 0.5 mhlikelihood
    // in R
    // log10(prod(0.5*(e/3) + 0.5*(1-e)))
    // -3.01156717
    MathTestUtils.assertAlmostEqual(mhlikelihood, -3.01156717)
  }
  test("Likelihood of simple model, all errors, no reference or mutant sites") {
    val m0likelihood = M0Model.logLikelihood("T", "G", all_c, None)
    val mflikelihood = MfmModel.logLikelihood("T", "G", all_c, None)
    val mhlikelihood = MHModel.logLikelihood("T", "G", all_c, None)

    // Assuming f = 0 or 0.5, mflikelihood or m0likelihood
    // in R:
    // options(digits=10)
    // p = seq(30,39)
    // e = 10^(-p/10)
    // log10(prod(e/3))
    // > -39.27121
    MathTestUtils.assertAlmostEqual(m0likelihood, -39.27121255)
    MathTestUtils.assertAlmostEqual(mflikelihood, m0likelihood)
    MathTestUtils.assertAlmostEqual(mhlikelihood, m0likelihood)

  }

  test("Likelihood of simple model, all mutant, no error or reference sites") {
    val m0likelihood = M0Model.logLikelihood("A", "C", all_c, None)
    // Assuming f = 0, mflikelihood or m0likelihood
    // in R:
    // options(digits=10)
    // p = seq(30,39)
    // e = 10^(-p/10)
    // log10(prod(e/3))
    // > -39.27121
    MathTestUtils.assertAlmostEqual(m0likelihood, -39.27121255)

    val mhlikelihood = MHModel.logLikelihood("C", "A", all_c, None)
    // assuming f = 0.5 mhlikelihood
    // in R
    // log10(prod(0.5*(1-e) + 0.5*(e/3)))
    // -3.01156717
    MathTestUtils.assertAlmostEqual(mhlikelihood, -3.01156717)

  }

  test("Likelihood/log10 likelihood of small mutant, first 3 reads (30,31,32) MfModel") {
    val mflikelihood = MfmModel.logLikelihood("C", "A", some_muts, None)
    val m03likelihood = MfmModel.logLikelihood("C", "A", some_muts, Some(0.3))
    val mhlikelihood = MHModel.logLikelihood("C", "A", some_muts, None)
    val m0likelihood = M0Model.logLikelihood("C", "A", some_muts, None)
    // f should be 0.3
    // in R
    // pm = seq(30,32)
    // pr = seq(33,39)
    // em = 10^(-pm/10)
    // er = 10^(-pr/10)
    // log10(prod(0.3*(1-em) + 0.7*(em/3))*prod( 0.3*(er/3) + 0.7*(1-er)))
    // -2.653910268
    MathTestUtils.assertAlmostEqual(mflikelihood, -2.653910268)
    MathTestUtils.assertAlmostEqual(m03likelihood, mflikelihood) //same if f properly calculated here

    // the no mutants model in R
    // log10(prod(0*(1-em) + 1*(em/3))*prod( 0*(er/3) + 1*(1-er)))
    // -10.73221105
    MathTestUtils.assertAlmostEqual(m0likelihood, -10.73221105)

    // the heterozygous site model in R:
    // log10(prod(0.5*(1-em) + 0.5*(em/3))*prod( 0.5*(er/3) + 0.5*(1-er)))
    // -3.01156717
    MathTestUtils.assertAlmostEqual(mhlikelihood, -3.01156717)
    val prob_mutant = MutectLogOdds.logOdds("C", "A", some_muts, None)
    MathTestUtils.assertAlmostEqual(prob_mutant, mflikelihood - m0likelihood)
  }

}
