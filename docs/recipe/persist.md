# Persistence Time Analysis

Computes the persistence time, from the continuous time correlation function,
following the definition by Luzar[@2000_Luzar] for hydrogen bonds.

## Usage

```bash
tame persist traj.dump -t '3,4:SSP:3.0,4.0'
```

## Options

| Option [shorthand] | Default | Description                                      |
|--------------------|---------|--------------------------------------------------|
| `-top`             | `None`  | Topology file                                    |
| `-dt`              | `1`     | Time step of the trajectory [ps]                 |
| `--seg [-s]`       | `5000`  | Segment length [ps]                              |
| `--format [-f]`    | `auto`  | Trajectory format                                |
| `--max-dt`         | `20`    | Time correlation function window                 |
| `--tags [-t]`       | `'1'`   | Tags definitions, see [below](##Tag definitions) |
| `--tcf-out`        | `'tcf'` | Time correlation function output                 |

## Tag definitions

Several definitions exist for the persistence time:

- The Impey, Madden and McDonald (IMM) definition uses a tolerance time to
  ignore transient escapes.[@1983_ImpeyMaddenEtAl]
- The static state picture (SSP) instead uses two cutoffs to define the "stable
  reactants" and "stable products".[@2008_LaageHynes]

\bibliography
