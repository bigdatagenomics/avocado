# avocado #

### Version parent ###
* ISSUE [95](https://github.com/bigdatagenomics/avocado/pull/95): Updating ADAM to 0.12.2-SNAPSHOT with migration to bdg-formats 0.2.0.
* ISSUE [93](https://github.com/bigdatagenomics/avocado/pull/93): SNAP streaming changes
* ISSUE [92](https://github.com/bigdatagenomics/avocado/pull/92): Adding a few fixes for running on a cluster.
* ISSUE [91](https://github.com/bigdatagenomics/avocado/pull/91): Adding fast aligner.
* ISSUE [90](https://github.com/bigdatagenomics/avocado/pull/90): Updating active region threshold.
* ISSUE [89](https://github.com/bigdatagenomics/avocado/pull/89): Adding release scripts.
* ISSUE [88](https://github.com/bigdatagenomics/avocado/pull/88): Removing duplicate write of header.
* ISSUE [87](https://github.com/bigdatagenomics/avocado/pull/87): Cleaning up avocado build.
* ISSUE [85](https://github.com/bigdatagenomics/avocado/pull/85): Fixing dependency issues due to injection of javax.servlet.
* ISSUE [84](https://github.com/bigdatagenomics/avocado/pull/84): Fixing read group code after #285 merged in adam upstream.
* ISSUE [83](https://github.com/bigdatagenomics/avocado/pull/83): Updating for move over to the bdg-formats repository.
* ISSUE [82](https://github.com/bigdatagenomics/avocado/pull/82): Fixed IO deadlock when running SNAP as a read frontend.
* ISSUE [81](https://github.com/bigdatagenomics/avocado/pull/81): Fixing issues with writing BAMs to standard input for external read call
* ISSUE [80](https://github.com/bigdatagenomics/avocado/pull/80): Changed license header to new BDGenomics license header.
* ISSUE [78](https://github.com/bigdatagenomics/avocado/pull/78): Fixing build.
* ISSUE [76](https://github.com/bigdatagenomics/avocado/pull/76): Adding generic stub for external variant callers.
* ISSUE [77](https://github.com/bigdatagenomics/avocado/pull/77): Refactoring assembler into a general haplotype based caller.
* ISSUE [75](https://github.com/bigdatagenomics/avocado/pull/75): Upgrade to ADAM 0.10.1-SNAPSHOT stream
* ISSUE [74](https://github.com/bigdatagenomics/avocado/pull/74): Refactoring transition matrix to allow for cleaner configuration
* ISSUE [72](https://github.com/bigdatagenomics/avocado/pull/72): Adding code to align reads with SNAP.
* ISSUE [70](https://github.com/bigdatagenomics/avocado/pull/70): Update config and debug
* ISSUE [64](https://github.com/bigdatagenomics/avocado/pull/64): Add flanks to gather
* ISSUE [63](https://github.com/bigdatagenomics/avocado/pull/63): Updated pair haplotype scoring method
* ISSUE [62](https://github.com/bigdatagenomics/avocado/pull/62): Change locuspredicate to unique read mapped predicate
* ISSUE [61](https://github.com/bigdatagenomics/avocado/pull/61): Fixes to HMM to add reference padding.
* ISSUE [60](https://github.com/bigdatagenomics/avocado/pull/60): Remove is<Type> methods from variantcall
* ISSUE [59](https://github.com/bigdatagenomics/avocado/pull/59): Get contig length from adamcontig
* ISSUE [52](https://github.com/bigdatagenomics/avocado/pull/52): Refactor HMMAligner and add testing for HMMAligner
* ISSUE [56](https://github.com/bigdatagenomics/avocado/pull/56): Parametrize flanking sequence length for k-mer graph creation
* ISSUE [55](https://github.com/bigdatagenomics/avocado/pull/55): Simplify haplotype queue creation
* ISSUE [46](https://github.com/bigdatagenomics/avocado/pull/46): Refactor KmerGraph creation
* ISSUE [53](https://github.com/bigdatagenomics/avocado/pull/53): Use function from ADAM to create reference sequence
* ISSUE [51](https://github.com/bigdatagenomics/avocado/pull/51): Update refid to use refname instead
* ISSUE [50](https://github.com/bigdatagenomics/avocado/pull/50): Parameterize HMM
* ISSUE [49](https://github.com/bigdatagenomics/avocado/pull/49): Fixes to HMM/Haplotype code to correct indel calling.
* ISSUE [48](https://github.com/bigdatagenomics/avocado/pull/48): Fix scheme issue with HDFS
* ISSUE [47](https://github.com/bigdatagenomics/avocado/pull/47): Move haplotype to separate classes
* ISSUE [45](https://github.com/bigdatagenomics/avocado/pull/45): Change snp caller interface to use rods exclusively
* ISSUE [44](https://github.com/bigdatagenomics/avocado/pull/44): Homozygous ref and alt were flipped in terms of writing/saving
* ISSUE [43](https://github.com/bigdatagenomics/avocado/pull/43): Upgrade adam to 0.9.1-SNAPSHOT and Spark to 0.9.1
* ISSUE [42](https://github.com/bigdatagenomics/avocado/pull/42): Use sequence load instead of only adam load for reference data
* ISSUE [41](https://github.com/bigdatagenomics/avocado/pull/41): Update Scalariform to be consistent with ADAM
* ISSUE [40](https://github.com/bigdatagenomics/avocado/pull/40): Move all files to proper package subdirectory
* ISSUE [38](https://github.com/bigdatagenomics/avocado/pull/38): Updated avocado to move to org.bdgenomics, and to move to ADAM 0.8.1-SNAPSHOT
* ISSUE [37](https://github.com/bigdatagenomics/avocado/pull/37): Splitting assembler into basic components (HMM, de Brujin graph, etc.)
* ISSUE [28](https://github.com/bigdatagenomics/avocado/pull/28): Fixes to assembler to get simple SNP calling working.
* ISSUE [34](https://github.com/bigdatagenomics/avocado/pull/34): Add input stage to sample configs
* ISSUE [32](https://github.com/bigdatagenomics/avocado/pull/32): Remove ReadModifier.scala, which appears to be dead code.
* ISSUE [30](https://github.com/bigdatagenomics/avocado/pull/30): Remove lingering references to VariantType, which was removed from ADAM
* ISSUE [29](https://github.com/bigdatagenomics/avocado/pull/29): Updated variant callers to remove variant type field, as it has been removed in ADAM
* ISSUE [21](https://github.com/bigdatagenomics/avocado/pull/21): Adding configurable input stage to pipeline.
* ISSUE [26](https://github.com/bigdatagenomics/avocado/pull/26): Updating to point to ADAM 0.7.2 snapshot release.
* ISSUE [24](https://github.com/bigdatagenomics/avocado/pull/24): Two changes: Update to build against ADAM 0.7.1-SNAPSHOT
* ISSUE [23](https://github.com/bigdatagenomics/avocado/pull/23): [WIP] Genotype likelihood calc was flipped
* ISSUE [22](https://github.com/bigdatagenomics/avocado/pull/22): Build changes for new schema
* ISSUE [20](https://github.com/bigdatagenomics/avocado/pull/20): Issue #19 Change references from ADAMFastaNucleotideContig to ADAMNucleo...
* ISSUE [17](https://github.com/bigdatagenomics/avocado/pull/17): Added a filter that filters genotypes that have been called with a high strand bias.
* ISSUE [18](https://github.com/bigdatagenomics/avocado/pull/18): Reupgraded to Spark 0.9.0. Added config for assembler.
* ISSUE [16](https://github.com/bigdatagenomics/avocado/pull/16): Reverts upgrade to spark 0.9.0
* ISSUE [15](https://github.com/bigdatagenomics/avocado/pull/15): Upgraded avocado to spark 0.9.0 and scala 2.10.3.
* ISSUE [14](https://github.com/bigdatagenomics/avocado/pull/14): Adding indel realigner to pipeline.
* ISSUE [13](https://github.com/bigdatagenomics/avocado/pull/13): Working pipeline. Added simple first cut configuration file.
* ISSUE [12](https://github.com/bigdatagenomics/avocado/pull/12): Adding optimized DefaultPartitionSet, better documentation, and tests.
* ISSUE [11](https://github.com/bigdatagenomics/avocado/pull/11): Adding support in for new FASTA contig class in ADAM.
* ISSUE [10](https://github.com/bigdatagenomics/avocado/pull/10): Cleanup 0.0.2
* ISSUE [8](https://github.com/bigdatagenomics/avocado/pull/8): General cleanup to remove cruft from merging code.
* ISSUE [7](https://github.com/bigdatagenomics/avocado/pull/7): Converting project over to maven.
* ISSUE [6](https://github.com/bigdatagenomics/avocado/pull/6): Updating with new variant/genotype implementations from Adam PR #13.
* ISSUE [4](https://github.com/bigdatagenomics/avocado/pull/4): Calling variants from local assembly
* ISSUE [3](https://github.com/bigdatagenomics/avocado/pull/3): EMforAlleles contains two main methods: emForMAFs and emForAFS.
* ISSUE [5](https://github.com/bigdatagenomics/avocado/pull/5): Cleaning up with new ADAMRod.
* ISSUE [2](https://github.com/bigdatagenomics/avocado/pull/2): Adding read filtering.
* ISSUE [1](https://github.com/bigdatagenomics/avocado/pull/1): Modified to allow VCF input containing MAFs for more advanced calling
