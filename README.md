avocado
=======

[![Coverage Status](https://coveralls.io/repos/github/bigdatagenomics/avocado/badge.svg?branch=master)](https://coveralls.io/github/bigdatagenomics/avocado?branch=master)

# A Variant Caller, Distributed

This README represents the TL;DR docs for avocado. More detailed documentation
is hosted at [Read the Docs](http://bdg-avocado.readthedocs.io/).

# Who/What/When/Where/Why avocado?

Avocado is a distributed variant caller built on top of the [ADAM format and
APIs](http://www.github.com/bigdatagenomics/adam) and [Apache
Spark](http://spark.apache.org/). Avocado is an open source project and is
released under the [Apache 2.0 license](https://github.com/bigdatagenomics/avocado/blob/master/LICENSE).

Avocado can be used for single sample germline variant calling, trio calling,
and joint variant calling. Avocado has >99% SNP calling accuracy, and >96%
INDEL calling accuracy when paired with ADAM's INDEL realignment pipeline.
When run on a single 32 core machine, Avocado can call variants on a 60x
coverage whole genome sequencing (WGS) dataset in approximately 7 hours. By
using Apache Spark to scale across multiple machines, Avocado can process the
same WGS dataset in approximately 15 minutes when using 1,024 cores.

# How avocado?

## Building Avocado

Avocado uses [Maven](http://maven.apache.org/) to build. To build avocado, cd
into the repository and run "mvn package".

## Avocado binaries

Nightly builds of Avocado are available from the [OSS Sonatype
repository](https://oss.sonatype.org/content/repositories/snapshots/org/bdgenomics/avocado/).
Additionally, we make a Docker image available from [Quay](https://quay.io/repository/ucsc_cgl/avocado?tag=latest&tab=tags).

# License

ADAM is released under the [Apache License, Version 2.0](LICENSE.txt).

# Citing Avocado

Avocado has been described in a PhD thesis. To cite this thesis, please cite:

```
@article{nothaft17,
  title={Scalable Systems and Algorithms for Genomic Variant Analysis},
  author={Nothaft, Frank Austin},
  school = {EECS Department, University of California, Berkeley},
  uear = {2017},
  month = {Dec},
  URL = {http://www2.eecs.berkeley.edu/Pubs/TechRpts/2017/EECS-2017-204.html},
  number = {UCB/EECS-2017-204}
}
```

A preprint describing Avocado should be released by the end of January 2018.