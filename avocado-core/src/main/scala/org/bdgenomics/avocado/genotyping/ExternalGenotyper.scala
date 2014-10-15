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

import org.seqdoop.hadoop_bam.{ SAMRecordWritable, VariantContextWritable }
import java.io.{ File, OutputStream }
import java.lang.Runnable
import java.nio.file.Files
import java.util.concurrent.{ ExecutorService, Executors }
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.SparkContext._
import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.AlignmentRecord
import org.bdgenomics.adam.converters.VariantContextConverter
import org.bdgenomics.adam.models.{ SAMFileHeaderWritable, VariantContext => ADAMVariantContext, ReferencePosition }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.GenomicRegionPartitioner
import org.bdgenomics.adam.rdd.read.AlignmentRecordContext._
import org.bdgenomics.avocado.models.{ Observation, ReadObservation }
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import htsjdk.variant.variantcontext.{ VariantContext => BroadVariantContext }
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.samtools.{ BAMStreamWriter, SAMFileHeader, SAMRecord }
import scala.collection.JavaConversions._
import scala.util.Sorting

object ExternalGenotyper extends GenotyperCompanion {

  val genotyperName: String = "ExternalGenotyper"

  protected def apply(stats: AvocadoConfigAndStats,
                      config: SubnodeConfiguration): Genotyper = {
    // get command and partition count
    val cmd = config.getString("command")
    val partitions = config.getInt("partitions")
    val debug = config.getBoolean("debug", false)

    new ExternalGenotyper(stats.contigLengths,
      partitions,
      cmd,
      debug)
  }
}

class ExternalGenotyper(contigLengths: Map[String, Long],
                        numPart: Int,
                        cmd: String,
                        debug: Boolean) extends Genotyper with Logging {

  val companion: GenotyperCompanion = ExternalGenotyper

  class ExternalWriter(records: Array[(ReferencePosition, SAMRecordWritable)],
                       header: SAMFileHeader,
                       stream: OutputStream) extends Runnable {

    def run() {
      val writer = new BAMStreamWriter(stream)

      // set header - header is written on first writeAlignment call
      writer.setHeader(header)

      var count = 0
      log.info("Starting to write records.")

      // cram data into standard in and flush
      records.foreach(r => {
        if (count % 100000 == 0) {
          stream.flush()
          log.info("Have written " + count + " reads.")
        }

        count += 1

        writer.writeHadoopAlignment(r._2)
      })

      // finish and close
      writer.close()
      stream.flush()
      log.info("Wrote a total of " + count + " reads")
    }
  }

  def callVariants(iter: Iterator[(ReferencePosition, SAMRecordWritable)],
                   header: SAMFileHeaderWritable): Iterator[VariantContextWritable] = {

    // sort records
    val reads = iter.toArray
    Sorting.stableSort(reads, (kv1: (ReferencePosition, SAMRecordWritable), kv2: (ReferencePosition, SAMRecordWritable)) => {
      kv1._1.compare(kv2._1) == -1
    })
    log.info("have " + reads.length + " reads")

    // we can get empty partitions if we:
    // - have contigs that do not have reads mapped to them
    // - don't have unmapped reads in our dataset
    if (reads.length == 0) {
      // can't write a bam file of length 0 ;)
      Iterator()
    } else {
      // get temp directory for vcf output
      val tempDir = Files.createTempDirectory("vcf_out")

      // build snap process
      val cmdUpdated = cmd.replaceAll("::VCFOUT::",
        tempDir.toAbsolutePath.resolve("calls.vcf").toString.replaceAll("\\\\", "\\\\\\\\"))
      log.info("running on host " + java.net.InetAddress.getLocalHost().getHostName() + ": " + cmdUpdated)
      val pb = new ProcessBuilder(cmdUpdated.split(" ").toList)

      // redirect error and get i/o streams
      pb.redirectError(ProcessBuilder.Redirect.INHERIT)

      // start process and get pipes
      val process = pb.start()
      val inp = process.getOutputStream()

      // get thread pool with two threads
      val pool: ExecutorService = Executors.newFixedThreadPool(2)

      // build java list of things to execute
      pool.submit(new ExternalWriter(reads, header.header, inp))

      // wait for process to finish
      val result = process.waitFor()
      if (result != 0) {
        log.error("Process " + cmdUpdated + " exited with " + result)
        throw new Exception("Process " + cmdUpdated + " exited with " + result)
      }

      // shut down pool
      pool.shutdown()

      log.info("process completed")

      // get vcf data - no index
      val vcfFile = new VCFFileReader(new File(tempDir.toAbsolutePath.toString + "/calls.vcf"), false)
      val iterator = vcfFile.iterator()
      var records = List[VariantContextWritable]()

      // loop and collect records
      while (iterator.hasNext()) {
        val record = iterator.next()

        // wrap variant context and append
        val vcw = new VariantContextWritable()
        vcw.set(record)
        records = vcw :: records
      }

      // need to close iterator - another code smell, but is required by samtools
      try {
        iterator.close()
      } catch {
        case ex: Exception =>
          log.warn("VCFReader exited with error " + ex)
      }

      // convert back to iterator and return
      records.toIterator
    }
  }

  def genotype(observations: RDD[Observation]): RDD[ADAMVariantContext] = {
    if (debug) {
      log.info("variant input DebugString:\n" + observations.toDebugString)
    }

    // get initial partition count
    val partitions = observations.partitions.length

    // convert records from adam to sam/bam
    val (reads, header) = observations.flatMap(r => r match {
      case rr: ReadObservation => Some(rr.read)
      case _                   => None
    }).adamConvertToSAM()

    // broadcast header
    val hdrBcast = reads.context.broadcast(SAMFileHeaderWritable(header))

    if (debug) {
      log.info("have " + reads.count + " reads")
    }

    // key reads by position and repartition
    val readsByPosition = reads.keyBy(r => ReferencePosition(r.get.getReferenceName.toString, r.get.getAlignmentStart))
      .partitionBy(new GenomicRegionPartitioner(numPart, contigLengths))

    if (debug) {
      log.info("have " + readsByPosition.count + " reads after partitioning")
    }

    // map partitions to external program
    val variants = readsByPosition.mapPartitions(r => callVariants(r, hdrBcast.value))

    // convert variants to adam format, coalesce, and return
    val converter = new VariantContextConverter()
    val result = variants.flatMap(vc => converter.convert(vc.get)).coalesce(partitions, true)
    if (debug) {
      log.info("variant output DebugString:\n" + result.toDebugString)
    }
    result
  }
}
