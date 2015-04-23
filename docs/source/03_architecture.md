# Architecture

avocado provides basic workflow management, which makes it simple to define a
repeatable analysis pipeline through the use of configuration files. The avocado
pipeline roughly supports the following workflow:

1. Any number of read pre-processing stages (e.g., BQSR, indel realignment,
PCR deduplications)
2. A variant discovery stage
3. A genotyping stage
4. Any number of variant processing/filtering stages (e.g., filter on depth or
strand bias, fit a VQSR-style model)

Configurations are processed via the
[Apache Commons Configuration](http://commons.apache.org/proper/commons-configuration/)
library. avocado uses the
[property list](http://commons.apache.org/proper/commons-configuration/userguide/howto_properties.html)
format for specifying hierarchical configurations. An odd side effect of our approach to
configurations is that a "configuration" must be specified for all stages, even if we
are not configuring that stage. For example, if we are going to use the indel realigner
preprocessing stage, but we aren't going to override any of the default settings, we must
provide an empty configuration object for that stage. We provide several example configurations
under the `avocado-sample-config/` directory.

This section provides an overview to avocado's pipeline architecture that is geared towards
people using avocado to call variants. If you are a developer who is interested in _extending_
avocado, you should look at the section on [avocado internals](#avocado Internals).

In addition to the methods described in this section, avocado contains experimental hooks for
running external tools (e.g., aligners, or other variant callers) distributed on Spark as part
of a variant calling pipeline. We have omitted these methods for now, as they are not yet stable
and remain work in progress.

## Read Pre-processing

avocado supports several read pre-processing stages, most of which come from the
[ADAM](https://www.github.com/bigdatagenomics/adam) project.

### Duplicate Marking

The duplicate marking stage identifies and filters reads which are possible PCR duplicates.
The algorithm used for identifying possible PCR duplicates is described in @massie13. We
have evaluated our algorithm against the equivalent algorithm in
[Picard](https://www.github.com/broadinstitute/picard); we are fully concordant, with the
exception of chimeric reads. Picard is unable to evaluate whether chimeric reads are PCR
duplicates, while our algorithm properly processes chimeric read pairs.

For configuration, the duplicate marker is named `markDuplicates`, and takes no options.

### Base Quality Score Recalibration

During the sequencing process, the sequencer may misestimate the base quality scores due
to global biases. Base quality score recalibration (BQSR) trains a simple statistical model
that identifies and corrects poor base quality score estimates. This model is constructed by
associating each read base with a covariate that represents an expected bias. For example,
we may expect biases to be correlated with the local sequence context, or the position of
the base in the read. An empirical quality score is estimated by applying either a Yates or
Bayesian model to the number of bases in each covariate that are an alignment mismatch and
the total number of bases in the covariate.

For configuration, BQSR is referred to via the stage name of `recalibrateBaseQualities`.
This stage supports one option:

* `snpTable`: If provided, loads a VCF file that describes variant sites that are masked
during BQSR. By default, no file is loaded.

### Indel Realignment

During global alignment, reads that span indels can be locally misaligned due to indel scoring
penalties. If enough reads that cover an indel are misaligned, the indel may not be called, if
using non-haplotype based variant discovery methods. Indel realignment identifies locations where
an indel may appear, and then locally realigns the reads that overlap this site. For more
details on the indel realignment algorithm, see @massie13 and @nothaft15.

> **Note:** Indel realignment is not recommended if using an assembly-based variant discovery
> approach. Running both together _will not_ impact accuracy, but it _will_ decrease performance.

For configuration, the indel realigner is referred to via the stage name of `realignIndels`.
Currently, the indel realigner in avocado does not support any options.

### Coalesce Reads

This preprocessing stage has no functional impact on the reads, but allows a user to control
the number of partitions that the reads are spread across. In turn, this impacts the parallelism
that avocado can achieve, and thus the performance of the pipeline.

The coalescing stage is referred to by the name `coalesceReads`. Users must pass the integer
argument `partitions`, which is the number of partitions to coalesce to. For guidance on the
number of partitions to choose, see the
[Spark tuning guide](http://spark.apache.org/docs/latest/tuning.html#level-of-parallelism).

## Variant Discovery

The variant discovery stage turns reads into evidence for variants. Currently, we have a single
approach for variant discovery that uses the read alignments. We are in the process of adding
a [local reassembly method](https://github.com/bigdatagenomics/avocado/pull/127) for detecting
variants.

### Discover From Reads

This method creates observations directly from the aligned reads. For each read, it traverses
from the alignment start to end. At each reference position, a single observation is generated.
As a result, all bases in an insertion are "compressed" down to a single site. For an _n_ base
deletion, _n_ observations will be generated. At the end, all observations from each read are
unioned with "observations" from the reference genome. Specifically, we take each base from
the provided reference assembly contigs and create an observation that tags the reference base
with its position.

### Discover By Reassembling Regions

This method takes all of the reads that align to a region of a reference contig and locally
reassembles the genome in that site. Our reasembly method differs from _most_ reassembly-based
approaches. In the conventional reassembly formulation, reads are put into an assembly graph.
Putative local haplotypes are then extracted from the reassembly graph. These haplotypes are
then scored by realigning the reads against the reassembly graphs. In our approach, we directly
emit allele observations from a de Brujin reassembly graph. This eliminates the need to generate
haplotypes and to realign reads, which both improves performance and allows for a less constrained
reassembly/realignment process. The approach we use to do this is described in more detail in
the [algorithms section](#local-reassembler).

The reassembler is referred to as `ReassemblyExplorer`. It takes a single parameter:

* `kmerLength`: The length of _k_-mers to use when building the de Brujin assembly graph. By
default, this is set to 20.

## Genotyping

Currently, the main genotyping model supported by avocado is a biallelic model that allows for
joint analysis. We are working on a somatic genotyping model, as well as a hereditary
genotyping model.

### Biallelic Genotyping Model

The biallelic genotyping model is based on the model used in the
[Samtools mpileup](https://samtools.github.io) engine [@li11]. In this genotyping model,
we assume independence between sites and apply a Bayesian model to the observed alleles
at each genomic locus. This model makes use of the mapping and base qualities to compute
the likelihood of each genotype state (where `g` equals the number of reference alleles
at each position). If the `emitGVCF` parameter is set to true, this model will emits variant
calls at all sites, including locations where a homozygous reference genotype is called.

We incorporate a [binomial distribution](en.wikipedia.org/wiki/Binomial_distribution)
as the prior distribution at each site, where the sample ploidy and major allele frequencies are
passed as parameters. We provide an EM algorithm to estimate the major allele frequency (MAF)
at each site from a group of samples; this EM algorithm is based on equations 13-17 from
@li11. Alternatively, if the user would like to simply use a population specific estimate
of the MAF, the EM algorithm can be disabled by setting `useEM` to `false`.

This variant calling model does not support multi-allelic sites. We plan to address this
with a more comprehensive, haplotype based variant calling method.

For configuration, the biallelic genotyper is named `BiallelicGenotyper`. It has the following
options:

* `ploidy`: The ploidy assumed at each site. Default value is 2.
* `useEM`: Whether or not to use an EM algorithm to find the MAF. Default is false.
* `emitGVCF`: Whether or not to emit confident reference genotype calls.
* `referenceFrequency`: The global expected frequency of reference alleles. If the EM algorithm
is used, this value is used to seed the EM loop. If the EM algorithm is disabled, this value is
used as the MAF in the prior distribution at each site.

There are three additional options specific to the EM algorithm. If the EM algorithm is being
used, the user must pick _at least_ one condition for terminating the EM loop. These include:

* `maxEMIterations`: The EM loop will stop running after _n_ iterations.
* `emSaturationThreshold`: The EM loop will stop running after a single step leads to a delta
smaller than this threshold value.

If both parameters are set, the EM algorithm will terminate once _either_ condition is achieved.

We have not designed the EM algorithm to ensure that the MAF will not saturate to either 0.0 or
1.0. Instead, we provide a final parameter to address the underflow problem:

* `emSaturationThreshold`: If the MAF estimate under/overflows, the site-specific MAF will clip
to either `emSaturationThreshold` or $1.0 - $`emSaturationThreshold`. By default, this parameter
is set to 0.001.

For additional detail about the math underlying the biallelic genotyper, refer to the [Biallelic
Genotyper] section.

### Simple Somatic Variant Calling

We currently have a stub for a simple somatic variant calling method. We plan to replace this
method soon with a reimplementation of MuTect [@cibulskis13]. The current method is not intended
for use; it is just a stub.

## Variant Post-processing

We provide a hook for post-processing variant calls after the variant calling/genotyping stage.
Most commonly, this will be used to filter out variant calls because they are likely spurious.
Right now, we support a limited set of post-processing stages, but we will be adding more after
a variant quality recalibration routine is added upstream in
[ADAM](https://www.github.com/bigdatagenomics/adam).

### Filter on Depth

This filter removes genotypes from sites with low read coverage. Users can specify an absolute
threshold, or a threshold that is relative to the average coverage across the sample.

For configuration, the depth filter is named `filterDepth`. It has the following two options:

* `absoluteDepth`: Filters out all genotypes without a depth greater than the absolute depth.
* `relativeDepth`: Filters out all genotypes with a depth lower than this ratio times the
average genome coverage.

At most, one of the `absoluteDepth` or `relativeDepth` options can be specified. If neither are
specified, we default to `relativeDepth = 0.75`.
