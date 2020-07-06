# Diffusion coefficient computations

There are generally two ways to calculate diffusion coefficients from
equilibrium molecular dynamics calculations: through the Green-Kubo relation or
through the corresponding Einstein relation:


Both methods are implemented in TAME, which are accessed through the `mdc` (mean
 displacement correlation) and `vcf` (velocity correlation function) subcommand
 of `tame diff`.

It's worth noting that the computation from diffusion coefficients, especially
from the MDC, is subject to the time range used in fitting/integration.
Therefore, both recipes save the correlations in `--corr-out` output file.

## Selection of diffusion coefficients

TAME can compute multiple (self/distinct) diffusion coefficient at once. The
diffusion coefficients are selected with the `--tags` flag. Self diffusion
coefficients are selected with the atom types, distinct diffusion coefficients
are selected with `type1,type2`. For instance, `tame diff mdc 1 2 1,2` computes
the self diffusion of type `1`, type `2` and the distinct diffusion coefficient
between `1` and `2`.

## Paired contribution to distinct diffusion coefficients

TAME can separate the pairing contribution with the following definition of ion
pairs: "pairs of atoms that stays within the cutoff radius for at least a give
lifetime $\tau$ at the time end when calculating the correlation function". This
type of distinct diffusion coefficient can be specified with the tag `1,2:C{rc}`.
Where `rc` is the cutoff radius.

$$
D^{d, paired}_{\alpha\beta} = \lim_{t \rightarrow \infty} \frac{1}{6t} \left [ \frac{2}{N} \sum_i^N  \sum_{j\ne i}^N \langle \Delta\mathbf{r}_{i,\alpha}(t) \cdot \Delta\mathbf{r}_{j,\beta}(t) \cdot \prod_{t'<\tau} B( R_{ij}(t')<R_{cut} )\rangle \right ]
$$

Alternatively, the paired contribution can be defined as following, where the
paired ions are defined as "pairs of atoms that stays within the cutoff radius
for the time window when the correlation function is calculated". This options
is specified with the `--loose` argument.

$$
D^\mathrm{d, paired}_{\alpha\beta} = \lim_{t \rightarrow \infty} \frac{1}{6t} \left [ \frac{2}{N}\sum_i^N  \sum_{j\ne i}^N \langle \Delta\mathbf{r}_{i,\alpha}(t) \cdot \Delta\mathbf{r}_{j,\beta}(t) \cdot \prod_{t'<t} B( R_{ij}(t')<R_{cut})\rangle \right ]
$$
