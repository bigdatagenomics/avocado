/*
 * Copyright (c) 2013. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.berkeley.cs.amplab.avocado.calls.pileup

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup,Base,ADAMGenotype,VariantType}
import edu.berkeley.cs.amplab.adam.models.{ADAMRod, ADAMVariantContext}
import edu.berkeley.cs.amplab.avocado.utils.Phred
import scala.math.pow
import scala.collection.mutable.MutableList
import scala.collection.JavaConversions._

/**
 * Class to call SNPs. Implements the SNP genotyping methods described in:
 *
 * Li, Heng. "A statistical framework for SNP calling, mutation discovery, association mapping and population
 * genetical parameter estimation from sequencing data." Bioinformatics 27.21 (2011): 2987-2993.
 *
 * for a single sample. As we are taking in a single sample, we do not calculate a minor allele frequency (MAF); rather,
 * we look up the MAF in a provided dbSNP vcf file. This provides MAF as the attribute MAF.
 * NOTE: In calculations we use MAJOR allele freqency for consistency with SAMTools publications! 
 * 
 * At the current point in time, we assume that we are running on a diploid organism.
 */
class PileupCallSNPVCFForMAF(sc: SparkContext, fileName: String) extends PileupCallSimpleSNP {
  
  override val callName = "SNPVCFForMAF"

  val maf = loadMaf(fileName, sc)
  val bcastMaf = sc.broadcast(maf)
  
  def loadMaf (fileName: String, sc: SparkContext): Map[Long,Double] = {
    val minMaf = 0.001
    val mafs = sc.newAPIHadoopFile[org.apache.hadoop.io.LongWritable, fi.tkk.ics.hadoop.bam.VariantContextWritable, fi.tkk.ics.hadoop.bam.VCFInputFormat](fileName)
      .map(rec =>(rec._2.get().getStart().toLong, rec._2.get().getAttributeAsDouble("GMAF", minMaf)))
      .collect()
      .toMap
    
    return mafs
  }

  /**
   * Compensates likelihoods. Likelihoods are compensated by MAF data.
   *
   * @param likelihoods List of likelihood values.
   * @param pileup List of pileups at this position.
   * @return Likelihoods after compensation.
   */
  override def compensate(likelihoods: List[Double], pileup: List[ADAMPileup]): List[Double] = {
    var compensatedLikelihoods = List[Double]()

    // get broadcast MAF values
    val mafs = bcastMaf.value

    // compensate likelihoods by allele frequency - get from vcf
    // NOTE: For consistency with SAMtools docs, we use MAJOR allele frequeny = 1 - maf
    val position = pileup.head.getPosition()
    val m = 1.0 - mafs(position)
    compensatedLikelihoods = likelihoods(2) * pow(m, 2.0) :: compensatedLikelihoods
    compensatedLikelihoods = likelihoods(1) * 2.0 * m * (1.0 - m) :: compensatedLikelihoods
    compensatedLikelihoods = likelihoods(0) * pow((1.0 - m), 2.0) :: compensatedLikelihoods
  
    compensatedLikelihoods
  }
}


