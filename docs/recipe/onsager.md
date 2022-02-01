# Onsager coefficients

Computing the onsager coefficients (RDF).

## Usage

```bash
tame onsager -t '3,3 3,4' traj.dump
```

## Options

| Option [shorthand] | Default | Description                      |
|--------------------|---------|----------------------------------|
| `-top`             | `None`  | Topology file                    |
| `-dt`              | `1`     | Time step of the trajectory [ps] |
| `--seg [-s]`       | `5000`  | Segment length [ps]              |
| `--format [-f]`    | `auto`  | Trajectory format                |
| `--tags [-t]`      | `'1,1'` | Space separated tags             |
| `--corr-out`       | `corr`  | Correlation function output      |

## Notes

The function outputs the volume-normalized displacement correlation function:

$$ \mathrm{corr}(t) = \frac{\langle \sum_i \Delta r_{i,\alpha}(t) \cdot \sum_j
\Delta r_{j,\beta}(t)\rangle}{\langle V \rangle} $$

with $\alpha$ and $\beta$ specified by the tags.

!!! warning "Partially implemented"

    Currently, the command does not perform fitting and does not output the Onsager coefficients
    directly.
