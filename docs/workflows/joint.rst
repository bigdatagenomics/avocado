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

In gVCF mode, The ``Jointer`` command runs by identifying sites where an
alternate allele was seen in at least one of the input files. We then use ADAM's
`region join <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/api/joins/>`__
functionality to square off and jointly genotype all of the gVCF inputs. These
three stages (discovering variants, extracting reference models, and joint
genotyping) can also be run as standalone steps. Running ``Jointer`` without the
``-from_gvcf`` flag will run the standalone joint genotyping step. The other two
standalone stages are discussed below.

Extracting Variants From A Set of gVCFs
---------------------------------------

To extract variants from a collection of gVCFs, run the ``Jointer`` command with
the ``-extract_only`` flag. This will scan over all genotype records in the input
files, and will identify the locations where at least one copy of the alternate
allele was seen. When run in this mode, instead of saving our output as a VCF
file, we save an Apache Parquet table of ADAM Variant records. This output can
be exported to VCF using `ADAM's transformVariants
command <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/cli/actions/#transformvariants>`__.

Extracting Reference Models from a gVCF
---------------------------------------

gVCF files represent all sites where reads were observed. If no alternate allele
was observed, the gVCF will contain a symbolic non-reference allele at the site,
represented as ``<NON_REF>``. This generic alternate allele is used to represent
the confidence that the reference was seen at this site. Even if we did not see
an alternate allele, we may not be able to confidently call the reference allele
if the site is poorly covered, poorly mapped, or covered by low quality bases.

To use gVCFs for joint genotyping, we need to extract the variants we are
interested in, and rewrite any generic non-reference records that overlap the
variant as the variant itself. To run only this phase, you can pass the
``-variants_to_extract`` parameter, along with a path to the variants you are
interested in extracting, encoded in any format ADAM can load as variants. The
output will be saved as Apache Parquet formatted using ADAM's Genotype schema.
This data can be exported to VCF using `ADAM's transformGenotypes
command <http://adam.readthedocs.io/en/adam-parent_2.11-0.23.0/cli/actions/#transformgenotypes>`__.
