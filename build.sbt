import AssemblyKeys._

assemblySettings

name := "Avocado"

version := "0.0.1"

organization := "edu.berkeley.cs.amplab"

scalaVersion := "2.9.3"

classpathTypes ++= Set ("orbit")

//classpathTypes ~= (_ + "orbit" + "bundle")

net.virtualvoid.sbt.graph.Plugin.graphSettings

excludeFilter := "ReadFilterOnComplexity.scala"

//exportJars := true
//seq(com.github.retronym.SbtOneJar.oneJarSettings: _*)
//mainClass in oneJar := Some("edu.berkeley.cs.amplab.avocado.Avocado")


libraryDependencies ++= Seq(
  "org.spark-project" % "spark-core_2.9.3" % "0.7.3",
  "org.scalatest" %% "scalatest" % "1.9.1" % "test",
  "edu.berkeley.cs.amplab.adam" % "adam-format" % "0.5.0-SNAPSHOT",
  "edu.berkeley.cs.amplab.adam" % "adam-commands" % "0.5.0-SNAPSHOT",
  "args4j" % "args4j" % "2.0.23",
  "commons-lang" % "commons-lang" % "2.6"
)

libraryDependencies ++= Seq(
  "org.eclipse.jetty.orbit" % "javax.servlet" % "2.5.0.v201103041518" artifacts (Artifact("javax.servlet", "jar", "jar")
  )
)

resolvers ++= Seq(
  "Typesafe" at "http://repo.typesafe.com/typesafe/releases",
  "Scala Tools Snapshots" at "http://scala-tools.org/repo-snapshots/",
  "Spray" at "http://repo.spray.cc",
  "Local" at "file:///Users/fnothaft/.m2/repository",
  "massie-maven" at "http://www.cs.berkeley.edu/~massie/maven/",
  "apache" at "https://repository.apache.org/content/repositories/releases"
)

mergeStrategy in assembly := {
 case m if m.toLowerCase.endsWith("manifest.mf") => MergeStrategy.discard
 case m if m.toLowerCase.matches("meta-inf.*\\.sf$") => MergeStrategy.discard
 case "META-INF/services/org.apache.hadoop.fs.FileSystem" =>
MergeStrategy.concat
 case "reference.conf" => MergeStrategy.concat
 case "log4j.properties" => MergeStrategy.concat
 case _ => MergeStrategy.first
}

assemblyCacheOutput in assembly := false