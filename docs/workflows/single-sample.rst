Single Sample Variant Calling
=============================

To call variants using Avocado, use the `BiallelicGenotyper` command. This
command discovers possibly variant sites from a collection of reads. The
discovered sites are then genotyped using a biallelic probabilistic model.
The genotyping model is based off of the biallelic model used by the original
`Samtools mpileup <https://academic.oup.com/bioinformatics/article/27/21/2987/217423>`__
variant caller, but adds additional components for modeling a site with two
minor alleles, as well as reads that do not match any known allele. The full
genotyping model is described in `Chapter 7 of this thesis <https://www2.eecs.berkeley.edu/Pubs/TechRpts/2017/EECS-2017-204.pdf>`__.

To run the `BiallelicGenotyper`, you must provide two parameters:

-  The path to the input file (in any file format for reads that `ADAM can load
   <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/api/adamContext/>`__)
-  The path where the output should be saved as Parquet-encoded Genotypes.

You can configure how Avocado discovers variants to genotype, the genotyping
phase, and the hard filters that Avocado applies to the called genotypes.

Reads to Evaluate
-----------------

By default, Avocado only calls variants in autosomal regions. Avocado does this
by inspecting the names of the contigs that the reads are mapped to. Avocado
assumes that the contig names start with ``chr`` prefixes. If your reference build
does not have these prefixes (i.e., chromosome 1 is ``1``, instead of ``chr1``), you
need to pass the ``-is_not_grc`` option. To enable calling non-autosomal regions,
you can pass:

-  ``-keep_mitochondrial_chromosome``: to call variants on the mitochondrial
   chromosome.
-  ``-autosomal_only``: to call variants on the sex chromosomes.

Additionally, we apply several read quality filters. To disable these filters,
you can pass:

-  ``-keep_duplicate_reads``: To disable discarding reads that have been marked as
   a PCR duplicate.
-  ``-keep_non_primary``: To disable discarding reads that have non-primary
   alignments.
-  ``-min_mapping_quality``: To change the default mapping quality below which
   reads are filtered out (default is phred 10).

Variant Discovery
-----------------

Avocado treats a site as possibly variant if enough alternate alleles are seen
with high enough quality. If you already know which variants you want to call,
you can skip variant discovery by passing the ``-variants_to_call`` flag, with a
path that ADAM can load as variants. During variant discovery, you can tune the
variants that are discovered with the following flags:

-  ``-min_phred_to_discover_variant``: The minimum Phred quality for a read
   containing a variant allele to be considered a confident observation of that
   allele.
-  ``-min_observations_to_discover_variant``: The minimum number of confident
   observations of a variant allele for us to choose to score a variant.

Genotyping
----------

By default, Avocado only scores the sites we have identified as possible
variants. However, Avocado can also emit genotype likelihoods for sites where
no alternate allele was seen. This output is not banded by quality and thus is
equivalent to a `BP RESOLUTION gVCF <https://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf>`__.
To run in this mode, pass the ``-score_all_sites`` option on the command line.

The only parameters that the genotyping engine consumes are around copy number.
Unless otherwise specified, Avocado assumes that it is running on a diploid
sample. To set the default copy number to another value, pass the ploidy with
the ``-ploidy`` option. Additionally, copy number variants can be passed with the
``-cnv`` flag. The copy number variants should be described in the GFF format that
is compatible with the `DECA <https://github.com/bigdatagenomics/deca>`__ and
`XHMM <http://atgu.mgh.harvard.edu/xhmm/index.shtml>`__ copy number variant
callers.

Variant Filtration
------------------

Avocado applies hard filters separately to SNPs and INDELs. The following hard
filters can be applied:

-  ``min`` and ``max``:
  -  ``depth``: The read depth covering the site.
  -  ``rms_mapping_quality``: The `RMS <https://en.wikipedia.org/wiki/Root_mean_square>`__
     quality of reads mapped at the site.
-  ``min`` and ``max``, separately for ``het`` and ``hom`` genotypes:
  -  ``allelic_fraction``: The fraction of reads that support the major vs.
     minor allele. E.g., if 25 reads support the major (reference) allele and
     75 reads support the minor (alternate) allele, the site has an allelic
     fraction of 0.75. For ``hom`` genotypes, only ``min`` is supported.
-  ``min`` only, separately for ``het`` and ``hom`` genotypes:
  -  ``quality_by_depth``: The genotype quality divided by the read depth. Sites
     with a low quality-by-depth may be highly covered (leading to high quality),
     but with low quality reads.

The default values are subject to change, but are described on the
``BiallelicGenotyper``'s command line when ``-help`` is printed.

Translating to VCF
------------------

The ``BiallelicGenotyper`` saves its output as `Apache Parquet <http://parquet.apache.org/>`__,
formatted using `ADAM's Genotype schema <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/architecture/schemas/>`__.
To translate this representation back to VCF, you can use `ADAM's
transformGenotypes command <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/cli/actions/#transformgenotypes>`__,
or Avocado's `Jointer <#joint>`__ command. We recommend using the ``Jointer``
command, which will calculate variant ``QUAL`` for all genotyped sites.
