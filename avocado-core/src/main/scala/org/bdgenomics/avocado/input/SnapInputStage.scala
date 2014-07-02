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
package org.bdgenomics.avocado.input

import java.io.{ OutputStream, InputStream, InputStreamReader }
import java.util.concurrent.{ Future, ExecutorService, Executors, Callable }
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.hadoop.io.{ LongWritable, Text }
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat
import org.apache.spark.SparkContext._
import org.apache.spark.{ Logging, SparkContext }
import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.{ ADAMRecord, ADAMNucleotideContigFragment }
import org.bdgenomics.adam.converters.SAMRecordConverter
import org.bdgenomics.adam.io.InterleavedFastqInputFormat
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.models.{ RecordGroupDictionary, SequenceDictionary }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import net.sf.samtools.{ SAMFileReader, SAMRecord, SAMReadGroupRecord, SAMRecordIterator }
import scala.collection.JavaConversions._

private[input] object SnapInputStage extends InputStage {

  val stageName = "SnapInputStage"

  /**
   * Sets up and loads data using this input stage.
   *
   * @param inputPath Path to input files.
   * @param config Configuration for this input stage.
   * @param stats Global stats/config data.
   * @return Returns an RDD of ADAM reads.
   */
  def apply(sc: SparkContext,
            inputPath: String,
            config: SubnodeConfiguration,
            reference: RDD[ADAMNucleotideContigFragment]): RDD[ADAMRecord] = {

    // build list of options
    // -M is included by default because ADAM assumes alignment match in cigar
    var cmdOpts: List[Option[String]] = List(Some("-M"),
      Some("-"),
      Some("-bam"),
      Some("-o"),
      Some("-"),
      Some("-pairedInterleavedFastq"),
      Some(config.getString("indexDirectory")),
      Some("paired"),
      Some(config.getString("snapPath", "snap")))

    val numMachines = config.getInt("numMachines")
    val coresPerMachine = config.getInt("coresPerMachine")

    /* Implement the following flags for SNAP:
     * 
     *   -d   maximum edit distance allowed per read or pair (default: 15)
     *   -n   number of seeds to use per read
     *   -sc  Seed coverage (i.e., readSize/seedSize).  Floating point.  Exclusive with -n.  (default: 0.000000)
     *   -h   maximum hits to consider per seed (default: 16000)
     *   -t   number of threads (default is one per core)
     *   -b   bind each thread to its processor (off by default)
     *   -e   compute error rate assuming wgsim-generated reads
     *   -P   disables cache prefetching in the genome; may be helpful for machines
     *        with small caches or lots of cores/cache
     *   -F   filter output (a=aligned only, s=single hit only, u=unaligned only)
     *   -I   ignore IDs that don't match in the paired-end aligner
     *   -E   misalign threshold (min distance from correct location to count as error)
     *   -Cxx must be followed by two + or - symbols saying whether to clip low-quality
     *        bases from front and back of read respectively; default: back only (-C-+)
     *   -M   indicates that CIGAR strings in the generated SAM file should use M (alignment
     *        match) rather than = and X (sequence (mis-)match)
     *   -G   specify a gap penalty to use when generating CIGAR strings
     *   --hp Indicates not to use huge pages (this may speed up index load and slow down alignment)
     *   -D   Specifies the extra search depth (the edit distance beyond the best hit that SNAP uses to compute MAPQ).  Default 2
     *   -rg  Specify the default read group if it is not specified in the input file
     */

    if (config.containsKey("maximumEditDistance")) {
      cmdOpts ::= Some("-d")
      cmdOpts ::= Some(config.getInt("maximumEditDistance").toString)
    }

    if (config.containsKey("seedsPerRead") && !config.containsKey("seedCoverage")) {
      cmdOpts ::= Some("-n")
      cmdOpts ::= Some(config.getInt("seedsPerRead").toString)
    } else if (!config.containsKey("seedsPerRead") && config.containsKey("seedCoverage")) {
      cmdOpts ::= Some("-sc")
      cmdOpts ::= Some(config.getDouble("seedCoverage").toString)
    } else if (config.containsKey("seedsPerRead") && config.containsKey("seedCoverage")) {
      throw new IllegalArgumentException("Cannot set both seed coverage and seeds per read in SNAP config.")
    }

    if (config.containsKey("unpairedMaxHitsPerRead")) {
      cmdOpts ::= Some("-h")
      cmdOpts ::= Some(config.getInt("unpairedMaxHitsPerRead").toString)
    }

    if (config.containsKey("pairedMaxHitsPerRead")) {
      cmdOpts ::= Some("-H")
      cmdOpts ::= Some(config.getInt("pairedMaxHitsPerRead").toString)
    }

    if (config.containsKey("snapThreads")) {
      cmdOpts ::= Some("-t")
      cmdOpts ::= Some(config.getInt("snapThreads").toString)
    }

    if (config.getBoolean("computeErrorRate", false)) {
      cmdOpts ::= Some("-e")
    }

    if (config.getBoolean("disablePrefetching", false)) {
      cmdOpts ::= Some("-P")
    }

    if (config.containsKey("filterOutput")) {
      cmdOpts ::= Some("-F")
      cmdOpts ::= Some(config.getString("filterOutput"))
    }

    if (config.getBoolean("ignoreID", false)) {
      cmdOpts ::= Some("-I")
    }

    if (config.containsKey("misalignThreshold")) {
      cmdOpts ::= Some("-E")
      cmdOpts ::= Some(config.getInt("misalignThreshold").toString)
    }

    if (config.containsKey("clipFrom")) {
      cmdOpts ::= Some("-C" + config.getString("clipFrom"))
    }

    if (config.containsKey("gapPenaltyForCigar")) {
      cmdOpts ::= Some("-G")
      cmdOpts ::= Some(config.getInt("gapPenaltyForCigar").toString)
    }

    if (config.getBoolean("dontUseBigPages", false)) {
      cmdOpts ::= Some("--hp")
    }

    if (config.containsKey("extraSearchDepth")) {
      cmdOpts ::= Some("-D")
      cmdOpts ::= Some(config.getInt("extraSearchDepth").toString)
    }

    if (config.containsKey("defaultReadGroup")) {
      cmdOpts ::= Some("-rg")
      cmdOpts ::= Some(config.getString("defaultReadGroup"))
    }

    if (config.containsKey("insertSizeRangeHigh")) {
      cmdOpts ::= Some("-s")
      cmdOpts ::= Some(config.getInt("insertSizeRangeLow", 0).toString)
      cmdOpts ::= Some(config.getInt("insertSizeRangeHigh").toString)
    }

    val cmd: List[String] = cmdOpts.flatMap((o: Option[String]) => o).reverse

    // set up snap runner
    val runner = new SnapRunner(cmd)

    // read fastq in from file
    val fastqIn: RDD[(Void, Text)] = sc.newAPIHadoopFile(inputPath,
      classOf[InterleavedFastqInputFormat],
      classOf[Void],
      classOf[Text])

    val fastqAsStrings: RDD[String] = fastqIn.map(kv => kv._2.toString)
      .coalesce(numMachines)

    // necessary data for conversion
    val seqDictBcast = fastqAsStrings.context.broadcast(reference.adamGetSequenceDictionary)

    // align reads using snap
    val reads = fastqAsStrings.mapPartitionsWithIndex(runner.mapReads(_, _, seqDictBcast.value))

    reads.coalesce(numMachines * coresPerMachine, true)
  }
}

/**
 * This is a class that sets up and runs one instance of the SNAP short-read aligner on a
 * cluster node during a mapPartitions phase. We open pipes to and from SNAP so that data
 * doesn't hit disk. Currently, we spit out SAMRecords; at a later point in time, we will
 * update SNAP to add ADAM as output.
 *
 * @param cmd Command line arguments to run, split up at whitespace.
 */
private[input] class SnapRunner(cmd: List[String]) extends Serializable with Logging {

  log.info("Running SNAP as: " + cmd.reduce(_ + " " + _))

  class SnapWriter(fastqStrings: Array[String],
                   stream: OutputStream) extends Callable[List[SAMRecord]] {
    def call(): List[SAMRecord] = {

      // cram data into standard in and flush
      fastqStrings.foreach(s => stream.write(s.getBytes))

      // flush and close stream
      stream.flush()
      stream.close()

      // this is a hack; require compatibility with Callable, which needs return
      List[SAMRecord]()
    }
  }

  class SnapReader(stream: InputStream)
      extends Callable[List[SAMRecord]] {

    def call(): List[SAMRecord] = {

      // push std out from snap into sam file reader
      val reader = new SAMFileReader(stream)

      // get iterator from sam file reader
      val iter: SAMRecordIterator = reader.iterator()

      var records = List[SAMRecord]()

      // loop through our poor iterator
      while (iter.hasNext()) {
        records ::= iter.next()
      }

      // close iterator, because that makes god damn sense
      iter.close()

      records
    }
  }

  /**
   * Method to map reads using SNAP. Given a set of strings which correspond to Fastq
   * input data. Returns an iterator of SAMRecords that have been wrapped to make them
   * serializable. Intended to be used in Spark as part of a mapPartitions phase.
   *
   * @param idx Partition index (used for logging).
   * @param fastqStrings An iterator containing strings corresponding to FASTQ records.
   * @param dict Sequence dictionary describing reference contigs.
   * @return An iterator containing aligned reads.
   */
  def mapReads(idx: Int, fastqStrings: Iterator[String], dict: SequenceDictionary): Iterator[ADAMRecord] = {

    // build snap process
    val pb = new ProcessBuilder(cmd)

    // redirect error and get i/o streams
    pb.redirectError(ProcessBuilder.Redirect.INHERIT)

    // start process and get pipes
    log.info("Starting SNAP on partition " + idx + "...")
    val process = pb.start()
    val inp = process.getOutputStream()
    val out = process.getInputStream()

    // get thread pool with two threads
    val pool: ExecutorService = Executors.newFixedThreadPool(2)

    // build java list of things to execute
    val exec: java.util.List[Callable[List[SAMRecord]]] = List(new SnapWriter(fastqStrings.toArray, inp),
      new SnapReader(out))

    // launch writer and reader in pool
    val futures: List[Future[List[SAMRecord]]] = pool.invokeAll(exec)

    // wait for process shutdown
    process.waitFor()
    log.info("SNAP is done running on partition " + idx + ".")

    // get read data
    val rec = futures.map(_.get())
    val records = rec(1) // select the records from the reader

    /**
     * Gets the name of the read group associated with a read.
     *
     * @param r Read to get read group from.
     * @return Option containing read group name if available.
     */
    def getReadGroupName(r: SAMRecord): Option[String] = Option(r.getReadGroup) match {
      case Some(rgr: SAMReadGroupRecord) => Some(rgr.getReadGroupId)
      case None                          => None
    }

    // collect read group names
    val readGroupNames = records.flatMap(getReadGroupName)
      .toSeq
      .distinct
    val readGroups = new RecordGroupDictionary(readGroupNames)

    // necessary data for conversion
    val recordConverter = new SAMRecordConverter

    // convert reads to ADAM and return
    val newRecs = records.map(recordConverter.convert(_, dict, readGroups))

    // shut down pool
    pool.shutdown()

    newRecs.toIterator
  }
}
