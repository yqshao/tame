# TAME Recipes

In TAME, recipes are scripts for specific properties. All recipes can be
accessed from the `tame` command line interface (CLI).

## General usage

In general, the TAME commands takes a (or multiple) trajectory file as input.
For instance, the below command computes the Onsager coefficients between
elements of type 1 and 2.

``` bash
tame onsager traj.dump --seg 5000 --dt 5 --tag '1;1' --tag '2;2' --tag '1;2'
```

A handy feature is that all TAME recipes allows for the split of the
trajectories into segments of equal lengths, and perform the analysis on them
separately. For instance, the above command will output several output files:
`mcd_seg1.dat`, `mcd_seg2.dat`, ... the variance thereof serves as an error
estimation.

## File format

The files formats are handled with the [chemfiles] library, which supports
analysis on a range of formats including, see the chemfiles
[documentation][formats] for detalis. The format can be supplied in CLI with the
`-f/--format` flag, e.g.:

```
tame onsager traj.dump -f LAMMPS -funit real
```

**Note**:

[chemfiles]: https://chemfiles.org
[formats]: https://chemfiles.org/chemfiles/latest/formats.html
[units]:

## List of commands

Below is a list of recipes implemented in TAME:

| Command   | Description                            |
| --------- | -------------------------------------- |
| `jacf`    | current autocorrelation function       |
| `mdc`     | mean displacement correlation          |
| `onsager` | Onsager coefficients                   |
| `pmsd`    | polarization mean squared displacement |
| `persist` | persistent time                        |
| `rdf`     | radial distribution function           |
| `sfac`    | structure factor                       |
