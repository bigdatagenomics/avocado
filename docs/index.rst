Avocado User Guide
==================

Introduction
============

Avocado is a variant caller built on top of `Apache
Spark <https://spark.apache.org>`__ to allow rapid variant calling on
cluster/cloud computing environments. Avocado is built on
`ADAM's <https://github.com/bigdatagenomics/adam>`__ APIs, and achieves variant
calling accuracy that is similar to state-of-the-art tools while being able
to drop variant calling latency to approximately 15 minutes when running on
a 1,024 core cluster.

Workflows Supported
===================

Avocado is run through the `avocado-submit` command line:

.. code:: bash

    ./bin/avocado-submit
   
::
   
   Using SPARK_SUBMIT=/usr/local/bin/spark-2.2.1-bin-hadoop2.7/bin/spark-submit

   Usage: avocado-submit [<spark-args> --] <avocado-args> [-version]

   Choose one of the following commands:

   biallelicGenotyper : Call variants under a biallelic model
   discover : Discover variants in reads
   jointer : Joint call and annotate variants.
   mergeDiscovered : Merge variants discovered from reads of multiple samples
   reassemble : Reassemble reads to canonicalize variants
   trioGenotyper : Call variants in a trio under a biallelic model

The `avocado-submit` script follows the same conventions as the `adam-submit`
command line, whose documentation can be found
`here <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/cli/overview/>`__.
As a result, just like ADAM, Avocado can be deployed on a local machine, on
`AWS <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/deploying/cgcloud/#running-adam-on-aws-ec2-using-cgcloud>`__,
an in-house cluster running `YARN <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/deploying/yarn/>`__
or `SLURM <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/deploying/slurm/>`__,
or using `Toil <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/deploying/toil/>`__.

Avocado supports several workflows:

-  Single sample germline variant calling: Avocado's
   BiallelicGenotyper runs on a single sample at a time, and can generate both
   variants-only (VCF) and all-sites (gVCF) output.
-  Joint variant calling: Avocado supports jointly calling
   variants from a collection of gVCF-styled inputs.

Avocado also contains code to reassemble variants, and a pedigree variant
caller. However, this code is experimental and is thus unsupported.

For genotyping, Avocado uses a probabilistic model that assumes that sites are
biallelic. This model is derived from the biallelic model used by the Mpileup
variant caller, but modified to better call multiallelic sites. When used with
the INDEL realigner from ADAM, Avocado has >99% accuracy when genotyping SNPs,
and >96% accuracy when genotyping INDELs.

.. toctree::
   :caption: Installation
   :maxdepth: 2

   installation/source

.. toctree::
   :caption: Workflows
   :maxdepth: 2

   workflows/single-sample
   workflows/joint

* :ref:`genindex`
* :ref:`search`
