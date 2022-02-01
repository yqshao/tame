# Radial ditribution function

Computing the radial distribution function (RDF).

## Usage

```bash
tame rdf -t '3,3 3,4' traj.dump
```

## Options

| Option [shorthand] | Default | Description                      |
|--------------------|---------|----------------------------------|
| `-top`             | `None`  | Topology file                    |
| `-dt`              | `1`     | Time step of the trajectory [ps] |
| `--seg [-s]`       | `5000`  | Segment length [ps]              |
| `--format [-f]`    | `auto`  | Trajectory format                |
| `--rmax`           | `30.0`  | Max radial distance              |
| `--rbin`           | `0.1`   | Radial distance bin size         |
| `--tag [-t]`       | `'1,1'` | Space separated tags             |
| `--rdf-out`        | `rdf`   | RDF output                       |
