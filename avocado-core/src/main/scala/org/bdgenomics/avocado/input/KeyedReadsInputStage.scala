package org.bdgenomics.avocado.input

import htsjdk.samtools.ValidationStringency
import org.apache.avro.Schema
import org.bdgenomics.adam.instrumentation.Timers._
import org.bdgenomics.formats.avro.{ AlignmentRecord, NucleotideContigFragment }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.ADAMContext
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.apache.hadoop.fs.{ FileStatus, FileSystem, Path }

private[input] object KeyedReadsInputStage extends InputStage {

  val stageName = "KeyedReads"

  /**
   * Sets up and loads data using this input stage.
   *
   * @param inputPath Path to input files.
   * @param config Configuration for this input stage.
   * @param reference RDD containing reference information.
   * @return Returns an RDD of ADAM reads.
   */
  def apply(sc: SparkContext,
            inputPath: String,
            config: SubnodeConfiguration,
            reference: RDD[NucleotideContigFragment]): RDD[AlignmentRecord] = {

    println("Loading reads in from " + inputPath)

    val path: Path = new Path(inputPath)
    var paths: Map[String, Path] = null
    if (inputPath.equals(config.getString("normalFileName"))) {
      paths = Map("normal" -> path)
    } else {
      paths = Map("tumor" -> path)
    }

    loadKeyedReads(sc, paths)
  }

  def loadKeyedReads(sc: SparkContext, paths: Map[String, Path]): RDD[AlignmentRecord] = {
    val alignments = paths map { case (k, p) => (k, sc.loadAlignments(p.toString).map(a => getAlignment(a, k))) }
    val alignedReads: Seq[RDD[AlignmentRecord]] = alignments.values.toList
    val union = sc.union(alignedReads)
    union
  }

  def getAlignment(a: AlignmentRecord, key: String) = {
    a.setRecordGroupDescription(key)
    a
  }

}
