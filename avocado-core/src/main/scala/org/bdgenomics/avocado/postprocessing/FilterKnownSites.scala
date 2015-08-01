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
package org.bdgenomics.avocado.postprocessing

import java.io.File
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{
  ReferencePosition,
  SnpTable,
  VariantContext
}
import org.bdgenomics.adam.rich.RichVariant
import org.bdgenomics.avocado.stats.AvocadoConfigAndStats
import org.bdgenomics.formats.avro.Genotype

object FilterKnownSites extends PostprocessingStage {

  val stageName = "filterKnownSites"

  def apply(rdd: RDD[VariantContext],
            stats: AvocadoConfigAndStats,
            config: SubnodeConfiguration): RDD[VariantContext] = {

    val snpTable = SnpTable(new File(config.getString("knownSiteFile")))
    val datasetName = if (config.containsKey("datasetName")) {
      Some(config.getString("datasetName"))
    } else {
      None
    }

    // by default, we _discard_ known sites
    val keepKnowns = config.getBoolean("keepKnowns", false)

    // build and apply filter
    val filter = new FilterKnownSites(datasetName, snpTable, keepKnowns)
    filter.filter(rdd)
  }
}

class FilterKnownSites(dataset: Option[String],
                       siteTable: SnpTable,
                       keepKnowns: Boolean) extends VariantAttributeFilter {

  val filterName = if (keepKnowns) {
    dataset.fold("KNOWN_SITE")(d => "SITE_IN_%s".format(d))
  } else {
    dataset.fold("UNKNOWN_SITE")(d => "SITE_NOT_IN_%s".format(d))
  }

  def filterFn(v: RichVariant): Boolean = {
    val in = siteTable.contains(ReferencePosition(v.variant.getContig.getContigName,
      v.variant.getStart))
    if (keepKnowns) {
      in
    } else {
      !in
    }
  }
}
