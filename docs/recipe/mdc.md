# Mean Displacement Correlation

This command computes the mean displacement correlation (MDC) for the
determination of diffusion coefficients.

When a single species is specified in the tag, the mean squared displacement
(MSD) of individual particles are computed which correspond to self diffusion
coefficients $D^{\mathrm{s}}$. When two are specified, the MDC between two
displacements are computed which correspond to distinct diffusion coefficients
$D^{\mathrm{d}}$.

$$
\begin{aligned}
D^\mathrm{s}_{\alpha} &= \lim_{t \to \infty}
\frac{1}{3tN} \left [ \sum_i \left\langle
  ||\Delta\mathbf{r}_{i,\alpha}(t)||^2
\right\rangle  \right]\\
D^\mathrm{d}_{\alpha\beta} &= \lim_{t \to \infty}
\frac{1}{3tN} \left [ \sum_i \sum_{j\ne i}  \left\langle
  \Delta\mathbf{r}_{i,\alpha}(t) \cdot \Delta\mathbf{r}_{j,\beta}(t)
\right\rangle \right]
\end{aligned}
$$

In addition, distinct diffusion can be combined with a definition of persistence
to get the pairing contribution of the diffusion coefficients according to the
following equation:[@2021_GudlaShaoEtAl]

$$ D^\mathrm{d, pairing}_{\alpha\beta} = \lim_{t \to \infty}
\frac{1}{3tN} \left [ \sum_i \sum_{j\ne i}  \left\langle
\Delta\mathbf{r}_{i,\alpha}(t) \cdot \Delta\mathbf{r}_{j,\beta}(t) \cdot
f(r_{ij};s))\right\rangle \right] $$

, where $f(r_{ij}; s)$ is a function related to the persistence time of the
pair, of which different definitions exist, to compute the distinct diffusion
coefficient for pairs, use a tag such as `'3;4;SSP:3.0,4.0'`, available
definitions of the persistence function can be found in the [`tame
persist`](persist.md#Persistence time definitions) documentation. This tag also
outputs the persistence time (as the count of alive pairs) and the pairing
distribution at each $t$ value, where $s=t$. Care should be taken when
interpreting the values, the point estimation of $\mathrm{MDC}(t)/3tN$ or its
fitted value is only meaningful if the persistence time converges **and** the
MDC is linear to time.

## Usage

```bash
tame mcd -t '3' -t '4;4' -t '3;4;SSP:3.0,4.0' traj.dump
```

## Options

| Option [shorthand] | Default | Description                 |
| ------------------ | ------- | --------------------------- |
| `-top`             | `None`  | Topology file               |
| `-dt`              | `1`     | Time step of the trajectory |
| `--seg [-s]`       | `5000`  | Segment length              |
| `--format [-f]`    | `auto`  | Trajectory format           |
| `--max-dt`         | `20`    | Time window for MDC fitting |
| `--tags [-t]`      | `'1'`   | Tags definitions            |
| `--mdc-out`        | `mdc`   | MDC output                  |

\bibliography
