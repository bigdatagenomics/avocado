Joint Variant Calling
=====================

Avocado's ``Jointer`` command supports joint variant calling from `gVCF-styled
data <https://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf>`__.
The ``Jointer`` command can also be used to export `Apache
Parquet <https://parquet.apache.org>`__ Genotype data to VCF, and to joint
genotype a collection of samples who all scored the same set of variants.
Our joint variant calling approach is is described in `Chapter 7 of this
thesis <https://www2.eecs.berkeley.edu/Pubs/TechRpts/2017/EECS-2017-204.pdf>`__.

To run the ``Jointer`` command, you must provide two parameters:

-  The path to all input files to joint genotyping (to load multiple files,
   use `Hadoop's glob syntax <https://hail.is/docs/stable/hadoop_glob_patterns.html>`__.
-  The path to save the output to, as a VCF file.

To save the VCF file as a single file (instead of sharded output), pass the
``-single`` flag.

If run on a single sample, ``Jointer`` will calculate variant statistics (VCF
INFO column attributes) and qualities only. If run on multiple samples, the
``Jointer`` command will update the called genotypes using a binomial prior
that is informed by the observed allele frequency of the variant across all
samples with confident calls. If the input data for the multiple samples is in
gVCF format, pass the ``-from_gvcf`` flag.
