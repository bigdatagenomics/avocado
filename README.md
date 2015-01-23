avocado
=======

A Variant Caller, Distributed

This README represents the TL;DR docs for avocado. More detailed documentation is under the `docs` folder.

# Who/What/When/Where/Why avocado?

avocado is a distributed variant caller built on top of the [ADAM format and pipeline](http://www.github.com/bigdatagenomics/adam) and [Apache Spark](http://spark.incubator.apache.org/). Currently, avocado is a research project, that is being built up to production level quality. avocado is open source and is covered by the Apache license.

There are several reasons that we are building avocado:

* As the price of genetic sequencing drops, we will need to be able to process greater and greater sums of genetic data. Ideally, we'll also want to process it more quickly. By using best-of-breed distributed systems design techniques, we hope to build a system that can scale to satisfy these domains.
* There is a dearth of well maintained open-source variant calling systems out there. We hope to build a system on which bioinformaticians can quickly and easily implement, test, and iterate on new algorithms, without needing to build their own infrastructure.
* 
Avocado is currently in its infancy, but we hope to have something interesting to show very soon!

# How avocado?

## Building avocado

avocado uses [Maven](http://maven.apache.org/) to build. To build avocado, cd into the repository and run "mvn package".
