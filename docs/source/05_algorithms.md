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

When scoring a site, we use the most frequently occurring alternate allele _across all samples_
as the alternate allele.

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

A "no-call" genotype will be emitted if no reads are observed at the site.
