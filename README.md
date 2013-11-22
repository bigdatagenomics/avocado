avocado
=======

A Variant Caller, Distributed

# Who/What/When/Where/Why avocado?

avocado is a distributed variant caller built on top of the [ADAM format and pipeline](http://www.github.com/massie/adam) and [Apache Spark](http://spark.incubator.apache.org/). Currently, avocado is a research project coming out of the [UC Berkeley AMPLab](https://amplab.cs.berkeley.edu/). avocado is open source and is covered by the Apache license.

There are several reasons that we are building avocado:

* As the price of genetic sequencing drops, we will need to be able to process greater and greater sums of genetic data. Ideally, we'll also want to process it more quickly. By using best-of-breed distributed systems design techniques, we hope to build a system that can scale to satisfy these domains.
* There is a dearth of well maintained open-source variant calling systems out there. We hope to build a system on which bioinformaticians can quickly and easily implement, test, and iterate on new algorithms, without needing to build their own infrastructure.
* Additionally, we've got a few cool ideas here at Berkeley about improving variant calling. Watch this space for more!

Avocado is currently in its infancy, but we hope to have something interesting to show very soon!

# How avocado?

## Building avocado

avocado uses [sbt](http://www.scala-sbt.org/) to build. To build avocado, cd into the repository and run "sbt compile".

Known issues with building:
* There is a current issue where the first time that the project is built, it fails on the first build before succeeding on the second build. This is because of an issue with a local dependency, that should hopefully be fixed soon.
* And... that local dependency is Adam. To install Adam locally, check out Adam and do "mvn package install" (for build instructions, see [the Adam readme](https://github.com/massie/adam)). We are picky though: currently, avocado depends on an [Adam pull request that hasn't merged in yet](https://github.com/massie/adam/pull/12). This pull request is available as the repository head at [chartl/adam](https://github.com/chartl/adam).

Standalone JARs are built using [sbt assembly](https://github.com/sbt/sbt-assembly). Once you have compiled, you can build the JAR by running "sbt assembly". Currently, the initial JAR takes about 10 minutes to build, but YMMV.