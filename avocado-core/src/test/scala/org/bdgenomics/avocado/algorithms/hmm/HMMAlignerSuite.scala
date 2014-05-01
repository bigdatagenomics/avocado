package org.bdgenomics.avocado.algorithms.hmm

import org.scalatest.FunSuite

class HMMAlignerSuite extends FunSuite {

  test("test perfect match alignment") {
    val hmm = new HMMAligner
    val refSequence = "TCGATCGA"
    val testSequence = "TCGATCGA"
    val qualities = "FFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === false)
    assert(alignment.alignmentStateSequence === "========")
    assert(alignment.alignedReference === "TCGATCGA")
    assert(alignment.alignedSequence === "TCGATCGA")

  }

  test("test single-snp alignment") {
    val hmm = new HMMAligner
    val refSequence = "TCGATCGA"
    val testSequence = "TCGTTCGA"
    val qualities = "FFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === true)
    assert(alignment.alignmentStateSequence === "===X====")
    assert(alignment.alignedReference === "TCGATCGA")
    assert(alignment.alignedSequence === "TCGTTCGA")
  }

  test("test multi-snp alignment") {
    val hmm = new HMMAligner
    val refSequence = "TCGATCGA"
    val testSequence = "TGGATGGA"
    val qualities = "FFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === true)
    assert(alignment.alignmentStateSequence === "=X===X==")
    assert(alignment.alignedReference === "TCGATCGA")
    assert(alignment.alignedSequence === "TGGATGGA")
  }

  test("test single insertion alignment") {
    val hmm = new HMMAligner
    val refSequence = "TCGATCGA"
    val testSequence = "TCGAATCGA"
    val qualities = "FFFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === true)
    assert(alignment.alignmentStateSequence === "===I=====")
    assert(alignment.alignedReference === "TCG_ATCGA")
    assert(alignment.alignedSequence === "TCGAATCGA")
  }

  test("test homopolymer single insertion alignment") {
    val hmm = new HMMAligner
    val refSequence = "TACCAATGTAA"
    val testSequence = "TACCCAATGTAA"
    val qualities = "FFFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === true)
    assert(alignment.alignmentStateSequence === "==I=========")
    assert(alignment.alignedReference === "TA_CCAATGTAA")
    assert(alignment.alignedSequence === "TACCCAATGTAA")
  }

  test("test long insertion alignment") {
    val hmm = new HMMAligner
    val refSequence = "TCGATCGA"
    val testSequence = "TCGACCCCTCGA"
    val qualities = "FFFFFFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === true)
    assert(alignment.alignmentStateSequence === "====IIII====")
    assert(alignment.alignedReference === "TCGA____TCGA")
    assert(alignment.alignedSequence === "TCGACCCCTCGA")
  }

  test("test single deletion alignment") {
    val hmm = new HMMAligner
    val refSequence = "TCGATCGA"
    val testSequence = "TCGATGA"
    val qualities = "FFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === true)
    assert(alignment.alignmentStateSequence === "=====D==")
    assert(alignment.alignedReference === "TCGATCGA")
    assert(alignment.alignedSequence === "TCGAT_GA")
  }

  test("test mixed insertion and deletion alignment") {
    val hmm = new HMMAligner
    val refSequence = "TCGATCGA"
    val testSequence = "TCGAGGTCG"
    val qualities = "FFFFFFFFFF"

    val alignment = hmm.alignSequences(refSequence, testSequence, qualities)

    assert(alignment.hasVariants === true)
    assert(alignment.alignmentStateSequence === "====II===D")
    assert(alignment.alignedReference === "TCGA__TCGA")
    assert(alignment.alignedSequence === "TCGAGGTCG_")
  }

}
