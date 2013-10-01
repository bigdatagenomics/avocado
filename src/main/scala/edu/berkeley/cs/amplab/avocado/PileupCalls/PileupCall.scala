package edu.berkeley.cs.amplab.avocado

import parquet.filter.{RecordFilter, UnboundRecordFilter}
import parquet.column.ColumnReader
import parquet.filter.ColumnRecordFilter._
import parquet.filter.ColumnPredicates._

class PileupMatchesReference extends UnboundRecordFilter 
{
  def bind(readers: Iterable[ColumnReader]): RecordFilter = 
  {
    not (column ("op", notEqualTo (MATCH)).bind(readers))
  }
}

object StackupCalls 
{
  
  def callVote ()
  {
  }
  
}
