/*
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.bdgenomics.avocado.postprocessing.mutect

import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{ ReferenceRegion, ReferenceMapping }

import scala.reflect.ClassTag

class Classified[T](val value: T, val classes: Set[String]) extends Serializable {

  def classify(cls: String): Classified[T] = new Classified[T](value, classes + cls)
  def hasClasses(cls: String*): Boolean = cls.forall(classes.contains)
  def hasSomeClasses(cls: String*): Boolean = cls.exists(classes.contains)

}

class ClassifiedReferenceMapping[T](tMapping: ReferenceMapping[T]) extends ReferenceMapping[Classified[T]] {
  override def getReferenceName(value: Classified[T]): String =
    tMapping.getReferenceName(value.value)
  override def getReferenceRegion(value: Classified[T]): ReferenceRegion =
    tMapping.getReferenceRegion(value.value)
}

class ClassifiedRDDFunctions[T](val rdd: RDD[Classified[T]])(implicit kt: ClassTag[T]) extends Serializable {

  def filterByClasses(cls: String*): RDD[Classified[T]] = rdd.filter(_.hasClasses(cls: _*))
  def filterBySomeClasses(cls: String*): RDD[Classified[T]] = rdd.filter(_.hasSomeClasses(cls: _*))
  def classify(cls: String): RDD[Classified[T]] = rdd.map(_.classify(cls))
  def values(): RDD[T] = rdd.map(_.value)
}

object ClassifiedContext {

  implicit def mappingToClassifiedMapping[T](mapping: ReferenceMapping[T]): ReferenceMapping[Classified[T]] =
    new ClassifiedReferenceMapping[T](mapping)

  implicit def classifiedRDDToClassifiedRDDFunctions[T](rdd: RDD[Classified[T]])(implicit kt: ClassTag[T]): ClassifiedRDDFunctions[T] =
    new ClassifiedRDDFunctions[T](rdd)

  implicit def toClassifiedRDD[T](rdd: RDD[T])(implicit kt: ClassTag[T]): RDD[Classified[T]] =
    rdd.map(t => new Classified[T](t, Set()))
}
