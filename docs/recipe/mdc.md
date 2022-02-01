# Mean Displacement Correlation 

Computes the mean displacement correlation (MDC) for the determination of self
or distinct diffusion coefficients. The options thereof are detailed below:

## Usage

```bash
tame mcd -t '3 4 3,4:SSP:3.0,4.0' traj.dump
```

## Options

| Option [shorthand] | Default | Description                                      |
|--------------------|---------|--------------------------------------------------|
| `-top`             | `None`  | Topology file                                    |
| `-dt`              | `1`     | Time step of the trajectory [ps]                 |
| `--seg [-s]`       | `5000`  | Segment length [ps]                              |
| `--format [-f]`    | `auto`  | Trajectory format                                |
| `--max-dt`         | `20`    | Time window for MDC fitting [ps]                 |
| `--tags [-t]`      | `'1'`   | Tags definitions, see [below](##Tag definitions) |
| `--fit-min`        | `5`     | Minimal time for msd fit [ps]                    |
| `--fit-max`        | `20`    | Maximal time for msd fit [ps]                    |
| `--mdc-out`        | `mdc`   | MDC output                                       |
| `--diff-out`       | `diff`  | diffusion coefficient output                     |

## Tags definitions

This command is capable of computing different diffusion coefficients in batch,
to specify multiple diffusion coefficients, delimit their definitions with
spaces in a string. The command supports the following types of tags:

- `1`: A single species specifies the self diffusion coefficient of the species
- `1,2`: Two species (can be the same) delimited by a comma specifies the
  distinct diffusion coefficient defined between them.
- `1,2:SSP:3,4`: Extending the distinct diffusion coefficient with a
  [persistence time definition](persist.md#Persistence definitions) separates the
  coefficient into the contribution of pairing and non-pairing contributions,
  see below.

## Pairing contributions

TAME can separate the pairing contribution of the diffusion coefficients
according to the following definition:[@2021_GudlaShaoEtAl]

$$ D^\mathrm{d, pairing}_{\alpha\beta} = \lim_{t \rightarrow \infty}
\frac{1}{3tN} \left [ \sum_i \sum_{j\ne i} \langle
\Delta\mathbf{r}_{i,\alpha}(t) \cdot \Delta\mathbf{r}_{j,\beta}(t) \cdot
f(r_{ij};s))\rangle \right] $$

, where $f(r_{ij}; s)$ is a function related to the persistence time of the pair,
of which different definitions exist, to compute the distinct diffusion
coefficient for pairs, use a tag such as `'3,4:SSP:3.0,4.0'`, available
definitions of the persistence function can be found in the [`tame
persist`](persist.md#Persistence time definitions) documentation.

Note that this tag outputs both the persistence time (as the count of alive
pairs) and the pairing distribution at each $t$ value, where $s=t$, which is
computationally more efficient in TAME. Care should be taken when interpreting
the values, the point estimation of $\frac{\mathrm{MDC}}{3tN}$ or its fitted
value is only meaningful if the persistence time converges **and** the
MDC is linear to time.

\bibliography
