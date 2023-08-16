# Onsager coefficients

This command outputs the displacement time correlation function whose time
derivatives correspond to the Onsager phenomenological coefficients
$\Omega_{\alpha\beta}$ between two species $\alpha$ and $\beta$ (specified by
the tags):

$$
\begin{aligned}
\mathrm{corr}_{\alpha\beta}(t)
&= \left\langle
    \sum_{ij} \Delta r_{i,\alpha}(t) \cdot \Delta
    r_{j,\beta}(t)
    \right\rangle \\
\Omega_{\alpha\beta}
&= \lim_{t\to\infty} \frac{\mathrm{corr}_{\alpha\beta}(t)}{6k_{\mathrm{B}}TVN_{\mathrm{A}}t}
\end{aligned}
$$

The time correlation function will be printed to the `corr-out` files.

## Usage

``` bash
tame onsager -t '3;3' -t '3;4' ...
```

## Options

| Option [shorthand] | Default | Description                          |
| ------------------ | ------- | ------------------------------------ |
| `-top`             | `None`  | Topology file                        |
| `-dt`              | `1`     | Time step of the trajectory          |
| `--seg [-s]`       | `5000`  | Segment length                       |
| `--max-dt`         | `20`    | Time window for correlation function |
| `--format [-f]`    | `auto`  | Trajectory format                    |
| `--rcom`           | `None`  | Remove centre of mass motion         |
| `--tag [-t]`       | `'1;1'` | Tags for pairs of atom types         |
| `--corr-out`       | `corr`  | Correlation function output          |
