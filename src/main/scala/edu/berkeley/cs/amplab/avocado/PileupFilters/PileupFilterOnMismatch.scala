package edu.berkeley.cs.amplab.avocado

import spark.SparkContext
import spark.SparkContext._
import parquet.filter.{RecordFilter, UnboundRecordFilter}
import parquet.column.ColumnReader
import parquet.filter.ColumnRecordFilter._
import parquet.filter.ColumnPredicates._

class PileupMatchesReference extends UnboundRecordFilter 
{
  def bind(readers: Iterable[ColumnReader]): RecordFilter = 
  {
    or (column ("op", equalTo (MISMATCH)),
	column ("op", equalTo (MISMATCH_REVERSE_STRAND))).bind(readers))
  }
}

class PileupHasMismatches extends UnboundRecordFilter
{
  def bind(readers: Iterable[ColumnReader]): RecordFilter =
  {
    column ("referenceId", in ())
  }
}

class PileupFilterOnMismatch 
{
  
}
