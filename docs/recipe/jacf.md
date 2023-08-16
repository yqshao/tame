# Current Auto-Correlation Function

Computing the ionic conductivity of from the Green-Kubo formula, using the
current auto-correlation function (JACF).

## Usage

```bash
tame jacf -c 3 -a 4 traj.dump
```

## Options

| Option [shorthand] | Default | Description                                        |
|--------------------|---------|----------------------------------------------------|
| `-top`             | `None`  | Topology file                                      |
| `-dt`              | `1`     | Time step of the trajectory [ps]                   |
| `--seg [-s]`       | `5000`  | Segment length [ps]                                |
| `--format [-f]`    | `auto`  | Trajectory format                                  |
| `--max-dt`         | `20`    | Time window for JACF                               |
| `--c-type [-c]`    | `3`     | Cation definitions, see [below](##Ion definitions) |
| `--a-type [-a]`    | `4`     | Anion definitions, see [below](##Ion definitions)  |
| `--temp [-T]`      | `293`   | Temperature [K]                                    |
| `--jacf-out`       | `jacf`  | JACF output                                        |
| `--cond-out`       | `cond`  | Conductivity output                                |
