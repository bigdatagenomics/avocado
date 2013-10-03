package edu.berkeley.cs.amplab.avocado.calls.pileup

import scala.math.pow

object Phred {

  /**
   * Static method for converting PHRED score to a probability.
   * Formula used is:
   * P = 10^(-Q/10)
   *
   * @param[in] phredScore Integer PHRED quality score.
   * @return Correctness likelihood.
   */
  def phredToProbability (phredScore: Int): Double = {
    return pow (10.0, (-phredScore / 10.0))
  }  
}
