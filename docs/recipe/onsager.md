# Onsager coefficients

Computing the onsager coefficients (RDF).

## Usage

```bash
tame onsager -t '3,3 3,4' traj.dump
```

## Options

| Option [shorthand] | Default | Description                               |
| ------------------ | ------- | ----------------------------------------- |
| `-top`             | `None`  | Topology file                             |
| `-dt`              | `1`     | Time step of the trajectory [ps]          |
| `--seg [-s]`       | `5000`  | Segment length [ps]                       |
| `--max-dt`         | `20`    | Time window for correlation function [ps] |
| `--format [-f]`    | `auto`  | Trajectory format                         |
| `--rcom`           | `None`  | Remove COM (see [notes](#Notes))          |
| `--tags [-t]`      | `'1,1'` | Space separated tags                      |
| `--corr-out`       | `corr`  | Correlation function output               |

## Notes

!!! warning "Partially implemented"

    Currently, the command does not perform fitting and does not output the Onsager coefficients
    directly.

### Correlation function

The function outputs the displacement correlation function:

$$
\mathrm{corr}(t) = \langle \sum_i \Delta r_{i,\alpha}(t) \cdot \sum_j \Delta
r_{j,\beta}(t)\rangle
$$

with $\alpha$ and $\beta$ specified by the tags.

### `--rcom` option

When specified, remove the center of mass motion of the entire system. A string
must be supplied with the mass of **all** species, e.g. `--rcom 1:1.008,3:6.941,6:12.010`.
