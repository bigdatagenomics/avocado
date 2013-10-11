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

package edu.berkeley.cs.amplab.avocado

import parquet.hadoop.{ParquetOutputFormat, ParquetInputFormat}
import spark.SparkContext
import spark.SparkContext._
import org.apache.hadoop.mapreduce.Job
import parquet.avro.{AvroParquetOutputFormat, AvroWriteSupport, AvroReadSupport}
import java.lang.Iterable
import com.google.common.io.Files
import java.io.File

object Avocado {

  def main (args: Array[String])
  {
    if (args.size != 1) {
      println ("Usage: avocado <config>")
    } else {
      val pipeline = new Avocado (args (0))
      
      pipeline.run ()
    }
  }
}

class Avocado (config: String) {
  
  val sc = new SparkContext ()
  

  def filterReads ()
  {
  }

  def processReads ()
  {
  }

  def filterPileups ()
  {
  }

  def callVariants ()
  {
  }

}
