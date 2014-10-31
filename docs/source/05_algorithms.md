# Algorithm Details

This section provides additional detail on algorithms implemented inside of avocado.
For details on the _use_ of these algorithms, refer to the [Architecture] section.

## Biallelic Genotyper

Mathematically, the biallelic genotyper is based on the Samtools mpileup engine [@li11].
We apply a Bayesian model where the prior distribution is a
[binomial distribution](en.wikipedia.org/wiki/Binomial_distribution) where the site ploidy
and the major allele frequency are used as priors. By default, we use an EM algorithm to
estimate the major allele frequency, however this can be overridden to default to a
genome-wide estimate of the reference allele frequency (0.999 in humans).

### Likelihood Function

The likelihood function used is:

$$
\mathcal{L}(g) = \frac{1}{m^{l + j}} \prod_{i = 1}^l (g * \epsilon_i + (m - g) * (1 - \epsilon_i))
\prod_{i = l + 1}^{l + j}(g * (1 - \epsilon_i) + (m - g) * \epsilon_i)
$$

Where the variables are defined as follows:

* $g$: genotype state, which is equal to the number of _reference alleles_ at this location
* $k$: The number of reads at this site.
* $i$: The number of reads at this site that match the reference allele.
* $j$: The number of reads at this site that match the alt allele.
* $\epsilon_i$: the error probability of the observation from read $i$, which is equal to
$1 - \sqrt{\text{mapQ} * \text{baseQ}}$

### EM Algorithm

The EM algorithm that we use estimates the site-specific MAF ($\psi$) from the per-sample
genotype likelihoods at a site. In each iteration $i$, we evaluate:

$$
\psi^{(i)} = \frac{1}{M} \sum_{s = 1}^{n} \frac{\sum_{g = 0}^{m_s} g P_s(g | \psi)
}{\sum_{g = 0}^{m_s} P_s(g | \psi)}
$$

Where:

* $\psi$ is the site-specific major (reference) allele frequency
* $i$ is the iteration
* $n$ is the number of samples at the site, and $s$ is the index of each sample
* $m_s$ is the ploidy of sample $s$ and $M = \sum_{s = 1}^{n} m_s$ is the site ploidy across
all samples
* $P_s(g | \psi)$ is the posterior probability of genotype state $g$ in sample $s$. The posterior
is $p(g | \psi) \mathcal{L}_s(g)$ where $\mathcal{L}_s(g)$ is the [genotype likelihood
function](#likelihood-function) for sample $s$ and $p(g | \psi)$ is the binomial prior.

### Genotyping

Once we have estimated the site specific MAF, $\psi$, we call the genotype state $\hat{g}$ as:

$$
\hat{g} = \argmax_g P(g | \psi) \mathcal{L}(g), g \in \{0, \dots , m\}
$$

The genotype quality is the phred scaled normalized probability of this genotype state:

$$
\text{Qual}_{\hat{g}} = -10 \log_{10} \frac{P(\hat{g} | \psi) \mathcal{L}(\hat{g})}{\sum_{g = 0}^{m} P(g | \psi) \mathcal{L}(g)}
$$

When scoring a site, we use the most frequently occurring alternate allele as the alternate
allele. A "no-call" genotype will be emitted if no reads are observed at the site.

## Local Reassembler

The local assembler in `avocado` is a de Brujin graph based assembler. The assembler creates
a "reference threaded" assembly. In this approach, we label _k_-mers that appear in the
reference genome with their reference position. The advantage of this approach is that we can
determine whether a _k_-mer in the graph can be mapped to the reference genome or if it
represents a divergence from the reference (an alternate allele). Additionally, for bubbles
off of reference arcs, we can identify the exact allele represented by the bubble without
needing to enumerate haplotypes which we then align against the reference sequence. This
approach has the following benefits:

* For a region that contains $n$ variants, we only need to evaluate $n$ variant arcs, as
opposed to $n^2$ variant haplotypes.
* We can emit statistics for computing genotype likelihoods directly from the de Brujin
graph. This is more efficient than emitting all haplotypes, scoring haplotype pairs, and
then realigning reads to the top scoring haplotype pair.

### Formulation

First, we start with the traditional formulation of a de Brujin graph for sequence assembly:

* Each $k$-mer $s$ represents a $k$-length string, with a $k - 1$ length prefix given by
$\text{prefix}(s)$ and a length 1 suffix given by $\text{suffix}(s)$.
* We place a directed edge ($\rightarrow$) from $k$-mer $s_1$ to $k$-mer $s_2$ if
$\text{prefix}(s_1)^{\{1, k - 2\}} + \text{suffix}(s_1) = \text{prefix}(s_2)$.

Additionally, we add the following:

* There is a set $\mathcal{R}$ which contains all of the $k$-mers that are in the reference
genome that covers a certain region.
* If $k$-mer $s \in \mathcal{R}$, then the output of function $\text{refPos}(s)$ is defined.
This function provides us with the integer position of $s$ in the reference genome.
* For two $k$-mers $s_1, s_2 \in \mathcal{R}$, we can define a distance function
$\text{distance}(s_1, s_2) = | \text{refPos}(s_1) - \text{refPos}(s_2) |$.

As stated above, this condition assumes that there is a 1-dimensional reference (e.g., a single
assembled reference). In reality, genomic assemblies provide a 2-dimensional space, where the
second dimension is defined by the chromosome that we are on. In practice, we will reassemble
reads from a region which is restricted to a single chromosome.

To discover alternate alleles from this graph _without_ performing a search for all haplotypes,
we need to address the $k$-mers that are not in $\mathcal{R}$. These can be processed with
path-based methods. First, let us classify paths into two types:

* *Spurs:* A spur is a set $S$ of $n$ $k$-mers $\{s_1, \dots, s_n\}$ where _either_ $s_1$ or
$s_n \in \mathcal{R}$ and all other $k$-mers are $\not\in \mathcal{R}$, and where
$s_i \rightarrow s_{i + 1} \forall i \in \{1, \dots, n - 1\}$. If $s_1 \in \mathcal{R}$,
then $s_n$ must not have a successor. Alternatively, if $s_n \in \mathcal{R}$, than $s_1$ is
required to not have a predecessor.
* *Bubbles:* A bubble is a set $S$ of $n$ $k$-mers $\{s_1, \dots, s_n\}$ where both
$s_1$ and $s_n \in \mathcal{R}$ and all other $k$-mers are $\not\in \mathcal{R}$, and where
$s_i \rightarrow s_{i + 1} \forall i \in \{1, \dots, n - 1\}$.

Currently, we do not process spurs. Spurs typically result from sequencing errors near the
start or end of a read. Additionally, given a spur, we cannot put a constraint on what sort of
edit it may be from the reference.

However, it is easy to define constraints for bubbles so that we can determine what variant
has been seen at that site. For many variants, we can even _canonically_ determine the variant.
First, let us make three definitions:

* The _bubble length_ is the number of non-reference $k$-mers that appear in the bubble.
* The _gap length_ is the distance between $s_1$ and $s_n$ in the reference genome. Formally,
$l_{\text{gap}} = \text{distance}(s_1, s_n)$.
* The _allele length_ is the absolute difference between the bubble length and the gap length.

With these definitions, we can define the allelic composition of bubbles as follows:

* A bubble contains a single/multiple nucleotide variant (S/MNV) if $l_{\text{allele}} = 0$.
* A bubble contains a canonical insert if $l_{\text{bubble}} = l_{\text{gap}} + l_{\text{allele}}$.
* A bubble contains a canonical deletion if $l_{\text{bubble}} < k$, where $k$ is the $k$-mer size,
and $\text{distance}(s_1, s_n) > 0$.
* If a bubble does not meet any of the above constraints, it is a non-canonical indel.

We discuss how to process these alleles in [the next section](#reassembly-implementation). It
is also worth noting that this approach can be used for detecting structural variants. A simple
extension involves jointly processing multiple reference regions (e.g., $\mathcal{R}_1$ and
$\mathcal{R}_2$). If $s_1$ and $s_n$ are members of different reference sets, then this bubble
identifies a possible structural variant breakpoint.

### Reassembly Implementation

Reassembly is performed as a two step process. In the first step, we build a reference threaded
de Brujin graph using the sequence of reference region $\mathcal{R}$ and the reads from sample $n$.
Once we have built a reference threaded de Brujin graph, we then elaborate the graph and identify
the allelic content of all bubbles.

#### Graph Construction

A reference threaded de Brujin graph is constructed per sample. We start by extracting the
reference $k$-mers. This algorithm---\texttt{addReferenceKmers}---is tail recursive.

\begin{algorithm}
\caption{Incorporate Reference \emph{k}-mers: \texttt{addReferenceKmers(iter, pos, lastKmer, kmerMap)}}
\label{alg:build-reference-kmers}
\begin{algorithmic}
\STATE $iter \leftarrow$ iterator across $k$-mer strings from the reference
\STATE $pos \leftarrow$ current reference position
\STATE $lastKmer \leftarrow$ the last $k$-mer seen
\STATE $kmerMap \leftarrow$ a map from $k$-mer strings $\rightarrow$ objects
\IF{$iter \ne \emptyset$}
\STATE $ks \leftarrow iter$.next
\IF{$lastKmer \ne \emptyset$}
\STATE $newKmer \leftarrow $ \texttt{Kmer(} $ks$, $pos$, $predecessors = \{lastKmer\}$\texttt{)}
\STATE $lastKmer$.$successors \leftarrow lastKmer$.$successors + newKmer$
\ELSE
\STATE $newKmer \leftarrow $ \texttt{Kmer(} $ks$, $pos$\texttt{)}
\ENDIF
\STATE $newPos \leftarrow pos + 1$
\STATE $kmerMap \leftarrow kmerMap + (ks, newKmer)$
\STATE \texttt{addReferenceKmers(} $iter$, $newPos$, $newKmer$, $kmerMap$\texttt{)}
\ENDIF
\end{algorithmic}
\end{algorithm}

We use a similar algorithm for processing the $k$-mers from the read dataset. However, this
algorithm filters out $k$-mers that contain \texttt{N}s, checks to see whether the $k$-mer
sequence has been seen before, and also stores information about the mapping quality, base
quality, and read strand.

#### Graph Elaboration

To elaborate the graph (the `toObservations` class method of a `KmerGraph`), we rely on a
tail-recursive state machine that pushes state onto a stack at branch points. Our machine
has the following states:

* `R`, `Reference`: We are on a run of $k$-mers that are mapped to a position in $\mathcal{R}$.
* `A`, `Allele`: We are on a run of $k$-mers that have diverged off of the reference. We have a
divergence start point, but have not connected back to the reference yet. This could be either
a bubble or a spur. 
* `C`, `ClosedAllele`: We were on a run of $k$-mers that had diverged off of the reference, but
have just reconnected back to the reference and now know the start and end positions of the bubble,
as well as the alleleic content of the bubble.

The following state transitions are allowed:

* `R` $\rightarrow$ `R`: We are continuing along a reference run.
* `R` $\rightarrow$ `A`: We were at a reference mapped $k$-mer, and have seen a branch to a
non-reference $k$-mer.
* `A` $\rightarrow$ `A`: We are continuing along a non-reference run.
* `A` $\rightarrow$ `C`: We were on a non-reference run, and have just connected back to a
reference $k$-mer.
* `C` $\rightarrow$ `R`: We have just closed out an allele, and are back at a reference position.

We initialize the state machine to `R` with the first $k$-mer from $\mathcal{R}$. Per $k$-mer,
we evaluate the state transition per successor $k$-mer. If the successor set contains a single
$k$-mer, we continue to that state. If the successor set contains multiple $k$-mers, we choose
a successor state to transition to, and push the branch context of all other successor states
onto our stack. If the successor set is empty, we pop a branch context off of the stack, and
switch to that context. We stop once we reach a $k$-mer whose successor set is empty, and our
branch context stack is empty.

The implementation of the `R` and `A` states are trivial and largely amount to bookkeeping. The
`C` state though, is more complex. First, let us discuss how to canonicalize the canonical alleles:

* For a canonical S/MNV bubble, the SNV sites can be determined by calculating the Hamming distance
of the bubble sequence versus the reference sequence. Evidence for SNVs exist at any site
where an edit is made.
* For a canonical insert, sequence is inserted $k - 1$ bases from the bubble start. There will be
$l_{\text{allele}}$ bases inserted. We can derive the inserted bases by computing the bubble
sequence and trimming the first $k - 1$ bases.
* For a canonical deletion, we guarantee that the bubble sequence is derived from the start and
end of the reference sequences that the bubble overlaps. To canonicalize this, we can perform
a greedy match of the bubble sequence, where we match bases to the front of the deletion until
we see a mismatch, and then align all remaining bases to the end of the deletion.

For a complex variant, we take the bubble sequence and align it against the reference sequence
from the bubble gap using a [pairwise HMM](#hmm-based-realignment). From the HMM alignment, we
emit observations according to the observed alignment states. The HMM uses an infinite padding
penalty to ensure that the ends of the sequences are matched.

## HMM-based Realignment

We use a HMM-based approach for pairwise sequence alignment, as described in @durbin98. When
using the aligner (`org.bdgenomics.avocado.algorithms.hmm.HMMAligner`), it takes a transition
matrix that specifies the various penalties (match, mismatch, gap open, gap continue, padding)
to use when scoring. The transition matrix used should vary:

* For pairwise alignment of two _reads_, the penalties should be based on the error model of the
sequencing technology.
* For alignment of a read to a reference, or the pairwise alignment of reference sequences, the
penalties should be based on variant frequencies.

Our HMM aligner emits an alignment object. This object includes an alignment state representation,
which can be run length encoded.
