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
import org.apache.spark.storage.StorageLevel
import org.apache.spark.{ Logging, SparkContext }
import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.{ ADAMRecord, ADAMNucleotideContigFragment }
import org.bdgenomics.adam.converters.SAMRecordConverter
import org.bdgenomics.adam.io.InterleavedFastqInputFormat
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.models.{ RecordGroup, RecordGroupDictionary, SequenceDictionary }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import net.sf.samtools._
import scala.collection.JavaConversions._

private[input] object SnapInputStage extends InputStage with Logging {

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
    val debug = config.getBoolean("debug", false)

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
    val rg = if (config.containsKey("defaultReadGroup")) config.getString("defaultReadGroup") else "FASTQ"
    val runner = new SnapRunner(cmd, rg)

    // read fastq in from file
    val fastqIn: RDD[(Void, Text)] = sc.newAPIHadoopFile(inputPath,
      classOf[InterleavedFastqInputFormat],
      classOf[Void],
      classOf[Text])

    if (debug) {
      log.info("SNAP input DebugString:\n" + fastqIn.toDebugString)
    }
    val fastqAsText: RDD[Text] = fastqIn.coalesce(numMachines, false).map(_._2)
    if (debug) {
      log.info("SNAP coalesced input DebugString:\n" + fastqAsText.toDebugString)
    }

    // necessary data for conversion
    val seqDictBcast = fastqAsText.context.broadcast(reference.adamGetSequenceDictionary)

    // align reads using snap
    val reads = fastqAsText.mapPartitionsWithIndex(runner.mapReads(_, _, seqDictBcast.value))

    val result = reads.coalesce(numMachines * coresPerMachine * 4, true).persist(StorageLevel.DISK_ONLY_2)
    if (debug) {
      log.info("SNAP output DebugString:\n" + result.toDebugString)
    }
    log.info("produced " + result.count + " reads")
    result
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
private[input] class SnapRunner(cmd: List[String], rg: String)
    extends Serializable with Logging {

  log.info("Running SNAP as: " + cmd.reduce(_ + " " + _))

  class SnapWriter(fastqStrings: Iterator[Text],
                   stream: OutputStream) extends Callable[Iterator[SAMRecord]] {
    def call(): Iterator[SAMRecord] = {

      // cram data into standard in and flush
      log.info("writing FASTQ strings for SNAP")
      var n = 0
      fastqStrings.foreach(s => {
        stream.write(s.getBytes)
        n += 1;
        if (n % 100000 == 0) {
          stream.flush;
          log.info("wrote " + n + " FASTQ strings for SNAP")
        }
      })

      // flush and close stream
      stream.flush()
      stream.close()
      log.info("wrote all " + n + " FASTQ strings for SNAP")

      // this is a hack; require compatibility with Callable, which needs return
      Iterator[SAMRecord]()
    }
  }

  class SnapReader(stream: InputStream, process: Process, pool: ExecutorService)
      extends Callable[Iterator[SAMRecord]] {

    def call(): Iterator[SAMRecord] = {

      class WrapIterator(inner: Iterator[SAMRecord])
          extends Iterator[SAMRecord] {
        var n = 0
        def next(): SAMRecord = {
          n += 1
          if (n % 100000 == 0) {
            log.info("read " + n + " SAM records from SNAP")
          }
          var again = true
          var result = new SAMRecord(new SAMFileHeader)
          while (again) {
            again = false
            try {
              result = inner.next()
            } catch {
              case ex: Exception =>
                log.warn("exception reading file" + ex)
                again = true
            }
          }
          result
        }
        def hasNext(): Boolean = {
          val more = inner.hasNext()
          if (!more) {
            log.info("read all " + n + " SAM records from SNAP")
            /*
            val result = process.waitFor()
            if (result != 0) {
              log.error("SNAP exited with error code " + result)
            } else {
              log.info("SNAP process finished")
            }
            pool.shutdown()
            */
          }
          more
        }
      }

      // push std out from snap into sam file reader
      // todo: ensure it gets closed?
      new WrapIterator(new SAMFileReader(stream).iterator())
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
  def mapReads(idx: Int, fastqStrings: Iterator[Text], dict: SequenceDictionary): Iterator[ADAMRecord] = {

    // build snap process
    val pb = new ProcessBuilder(cmd)

    // redirect error and get i/o streams
    pb.redirectError(ProcessBuilder.Redirect.INHERIT)

    // start process and get pipes
    log.info("Starting SNAP on partition " + idx + " on host " + java.net.InetAddress.getLocalHost().getHostName())
    val process = pb.start()
    val inp = process.getOutputStream()
    val out = process.getInputStream()

    // get thread pool with two threads
    val pool: ExecutorService = Executors.newFixedThreadPool(1)

    // build collection of things to execute, fork into separate thread
    pool.submit(new SnapWriter(fastqStrings, inp))

    // get read data from SNAP stdout, shutdown process & pool at end
    val records = new SnapReader(out, process, pool).call()

    // necessary data for conversion
    val readGroupDict = new RecordGroupDictionary(List(new RecordGroup(rg, rg)))
    val recordConverter = new SAMRecordConverter

    // convert reads to ADAM and return
    log.info("converting to ADAM")
    val newRecs = records.map(recordConverter.convert(_, dict, readGroupDict))

    newRecs.toIterator
  }
}
