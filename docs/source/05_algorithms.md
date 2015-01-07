# Algorithm Details

This section provides additional detail on algorithms implemented inside of avocado.
For details on the _use_ of these algorithms, refer to the [Architecture] section.

## Biallelic Genotyper

Mathematically, the biallelic genotyper is based on the Samtools mpileup engine [@li11].
The likelihood function used is:

$$
\mathcal{L}(g) = \frac{1}{m^{i + j}} \prod_{i = 1}^l (g * \epsilon_i + (m - g) * (1 - \epsilon_i))
\prod_{i = l + 1}^{k - j}(g * (1 - \epsilon_i) + (m - g) * \epsilon_i)
$$

Where the variables are defined as follows:

* $g$: genotype state, which is equal to the number of _reference alleles_ at this location
* $k$: The number of reads at this site.
* $i$: The number of reads at this site that match the reference allele.
* $j$: The number of reads at this site that match the alt allele.
* $\epsilon_i$: the error probability of the observation from read $i$, which is equal to
$1 - \sqrt{\text{mapQ} * \text{baseQ}}$

When scoring a site, we use the most frequently occurring alternate allele as the alternate
allele.
