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
package org.bdgenomics.avocado.serialization

import com.esotericsoftware.kryo.Kryo
import org.apache.spark.serializer.KryoRegistrator
import org.bdgenomics.adam.serialization.ADAMKryoRegistrator

class AvocadoKryoRegistrator extends KryoRegistrator {

  private val akr = new ADAMKryoRegistrator()

  override def registerClasses(kryo: Kryo) {

    // register adam's requirements
    akr.registerClasses(kryo)

    // org.bdgenomics.avocado.util.TreeRegionJoin
    kryo.register(classOf[org.bdgenomics.avocado.util.Forest[org.bdgenomics.formats.avro.Variant]],
      new org.bdgenomics.avocado.util.ForestSerializer[org.bdgenomics.formats.avro.Variant, org.bdgenomics.adam.serialization.AvroSerializer[org.bdgenomics.formats.avro.Variant]](
        new org.bdgenomics.adam.serialization.AvroSerializer[org.bdgenomics.formats.avro.Variant]))
  }
}

