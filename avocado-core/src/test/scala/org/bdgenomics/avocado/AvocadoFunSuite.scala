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
package org.bdgenomics.avocado

import org.bdgenomics.utils.misc.SparkFunSuite
import java.net.URL
import java.nio.file.Files

trait AvocadoFunSuite extends SparkFunSuite {

  override val appName: String = "avocado"
  override val properties: Map[String, String] = Map(("spark.serializer", "org.apache.spark.serializer.KryoSerializer"),
    ("spark.kryo.registrator", "org.bdgenomics.adam.serialization.ADAMKryoRegistrator"),
    ("spark.kryoserializer.buffer", "4M"),
    ("spark.kryo.referenceTracking", "true"))

  /** FIXME: STOLEN FROM ADAMFunSuite, see ADAM/#1225 **/
  def resourceUrl(path: String): URL = ClassLoader.getSystemClassLoader.getResource(path)

  def tmpFile(path: String): String = Files.createTempDirectory("").toAbsolutePath.toString + "/" + path
}

