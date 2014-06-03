/*
 * Copyright (c) 2014. Regents of the University of California
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
package org.bdgenomics.avocado.input

import fi.tkk.ics.hadoop.bam.SAMRecordWritable
import java.io.{ OutputStream, InputStream, InputStreamReader }
import java.util.concurrent.{ Future, ExecutorService, Executors, Callable }
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.hadoop.io.{ LongWritable, Text }
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat
import org.apache.spark.SparkContext._
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.avro.{ ADAMRecord, ADAMNucleotideContigFragment }
import org.bdgenomics.adam.converters.SAMRecordConverter
import org.bdgenomics.adam.io.InterleavedFastqInputFormat
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.models.RecordGroupDictionary
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

    // for now, assume that fastq records are all 4 lines; jeremy is providing a better fix
    val fastqAsStrings: RDD[String] = fastqIn.map(kv => kv._2.toString)
      .coalesce(numMachines)

    // align reads using snap
    val alignedReads: RDD[SAMRecordWritable] = fastqAsStrings.mapPartitions(runner.mapReads)

    /**
     * Gets the name of the read group associated with a read.
     *
     * @param r Read to get read group from.
     * @return Option containing read group name if available.
     */
    def getReadGroupName(r: SAMRecordWritable): Option[String] = Option(r.get.getReadGroup) match {
      case Some(rgr: SAMReadGroupRecord) => Some(rgr.getReadGroupId)
      case None                          => None
    }

    // collect read group names
    val readGroupNames = alignedReads.flatMap(r => getReadGroupName(r))
      .distinct
      .collect
      .toSeq
    val readGroups = new RecordGroupDictionary(readGroupNames)

    // necessary data for conversion
    val recordConverter = new SAMRecordConverter
    val readGroupsBcast = alignedReads.context.broadcast(readGroups)
    val seqDictBcast = alignedReads.context.broadcast(reference.adamGetSequenceDictionary)

    // convert aligned reads into adam
    val readsAsAdam: RDD[ADAMRecord] = alignedReads.coalesce(numMachines * coresPerMachine, true)
      .map(r => recordConverter.convert(r.get, seqDictBcast.value, readGroupsBcast.value))

    readsAsAdam
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
private[input] class SnapRunner(cmd: List[String]) extends Serializable {

  println(cmd.reduce(_ + " " + _))

  class SnapWriter(fastqStrings: Iterator[String],
                   stream: OutputStream) extends Callable[List[SAMRecordWritable]] {
    def call(): List[SAMRecordWritable] = {

      var count = 0
      println("Starting to write records.")

      // cram data into standard in and flush
      fastqStrings.foreach(s => {
        if (count % 1000 == 0) {
          println("Have written " + count + " reads.")
          stream.flush()
        }

        count += 1

        stream.write(s.getBytes)
      })
      stream.flush()

      // this is a hack; require compatibility with Callable, which needs return
      List[SAMRecordWritable]()
    }
  }

  class SnapReader(stream: InputStream)
      extends Callable[List[SAMRecordWritable]] {

    def call(): List[SAMRecordWritable] = {
      println("Starting reader.")

      // push std out from snap into sam file reader
      val reader = new SAMFileReader(stream)

      // get iterator from sam file reader
      val iter: SAMRecordIterator = reader.iterator()

      var count = 0
      println("Starting to read records.")

      var writableRecords = List[SAMRecordWritable]()

      // loop through our poor iterator
      while (iter.hasNext()) {
        if (count % 1000 == 0) {
          println("Have read " + count + " aligned reads.")
        }

        count += 1
        val rw = new SAMRecordWritable()
        rw.set(iter.next())

        writableRecords ::= rw
      }

      // close iterator, because that makes god damn sense
      iter.close()

      writableRecords
    }
  }

  /**
   * Method to map reads using SNAP. Given a set of strings which correspond to Fastq
   * input data. Returns an iterator of SAMRecords that have been wrapped to make them
   * serializable. Intended to be used in Spark as part of a mapPartitions phase.
   *
   * @param fastqStrings Iterator containing strings corresponding to FASTQ records.
   * @return An iterator containing aligned reads wrapped as SAMRecordWritable.
   */
  def mapReads(fastqStrings: Iterator[String]): Iterator[SAMRecordWritable] = {
    // build snap process
    val pb = new ProcessBuilder(cmd)

    // redirect error and get i/o streams
    pb.redirectError(ProcessBuilder.Redirect.INHERIT)

    // start process and get pipes
    println("Starting SNAP...")
    val process = pb.start()
    println("SNAP is started.")
    val inp = process.getOutputStream()
    val out = process.getInputStream()
    println("Have pipes.")

    // get thread pool with two threads
    val pool: ExecutorService = Executors.newFixedThreadPool(2)

    // build java list of things to execute
    println("Starting reader and writer.")
    val exec: java.util.List[Callable[List[SAMRecordWritable]]] = List(new SnapWriter(fastqStrings, inp),
      new SnapReader(out))

    // launch writer and reader in pool
    val futures: List[Future[List[SAMRecordWritable]]] = pool.invokeAll(exec)

    // get read data
    val records = futures.map(_.get())
      .reduce(_ ++ _)
      .toIterator

    // shut down pool
    pool.shutdown()

    records
  }
}
