# Polarization Mean Squared Displacement

Computing the ionic conductivity of from the Green-Kubo formula in the Einstein
form, using the polarization mean squared displacement (PMSD).

## Usage

```bash
tame pmsd -c 3 -a 4 traj.dump
```

## Options

| Option [shorthand] | Default | Description                                        |
|--------------------|---------|----------------------------------------------------|
| `-top`             | `None`  | Topology file                                      |
| `-dt`              | `1`     | Time step of the trajectory [ps]                   |
| `--seg [-s]`       | `5000`  | Segment length [ps]                                |
| `--format [-f]`    | `auto`  | Trajectory format                                  |
| `--max-dt`         | `20`    | Time window for PMSD                               |
| `--c-type [-c]`    | `3`     | Cation definitions, see [below](##Ion definitions) |
| `--a-type [-a]`    | `4`     | Anion definitions, see [below](##Ion definitions)  |
| `--temp [-T]`      | `293`   | Temperature [K]                                    |
| `--fit-min`        | `5`     | Minimal time for msd fit [ps]                      |
| `--fit-max`        | `20`    | Maximal time for msd fit [ps]                      |
| `--pmsd-out`       | `pmsd`  | PMSD output                                        |
| `--cond-out`       | `cond`  | conductivity output                                |
