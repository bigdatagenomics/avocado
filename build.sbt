name := "Avocado"

version := "0.0.1"

organization := "edu.berkeley.cs.amplab"

scalaVersion := "2.9.3"

libraryDependencies ++= Seq(
  "org.spark-project" % "spark-core_2.9.3" % "0.7.3",
  "org.scalatest" %% "scalatest" % "1.9.1" % "test"
  "org.streum" %% "configrity-core" % "1.0.0"
)

resolvers ++= Seq(
  "Typesafe" at "http://repo.typesafe.com/typesafe/releases",
  "Scala Tools Snapshots" at "http://scala-tools.org/repo-snapshots/",
  "Spray" at "http://repo.spray.cc"
)
