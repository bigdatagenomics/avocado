# Installation

## Building avocado from Source

You will need to have [Maven](http://maven.apache.org/) installed in order to build avocado.

> **Note:** The default configuration is for Hadoop 2.2.0. However, avocado does not package
> either Spark or Hadoop, and should not need to be rebuilt to run on different Spark (after 1.0.0)
> or Hadoop versions.

```
$ git clone https://github.com/bigdatagenomics/avocado.git
$ cd avocado
$ export "MAVEN_OPTS=-Xmx512m -XX:MaxPermSize=128m"
$ mvn clean package
```

## Running avocado

avocado is packaged via [appassembler](http://mojo.codehaus.org/appassembler/appassembler-maven-plugin/)
and packages all necessary dependencies _except for_ Spark and Hadoop, which must be present on the
system (see below).

You might want to add the following to your `.bashrc` to make running avocado easier:

```
alias avocado-submit="${AVOCADO_HOME}/bin/avocado-submit"
```

`$AVOCADO_HOME` should be the path to where you have checked avocado out on your local filesystem. 
This alias wraps the `spark-submit` command to run avocado. You'll need
to have the Spark binaries on your system; prebuilt binaries can be downloaded from the
[Spark website](http://spark.apache.org/downloads.html). If you do not have a cluster, we recommend using
[Spark 1.1, and Hadoop 2.3 (CDH5)](http://d3kbcqa49mib13.cloudfront.net/spark-1.1.0-bin-hadoop2.3.tgz).

