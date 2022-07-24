# TAME Recipes 

In TAME, recipes are scripts for specific properties. All recipes can be
accessed from the `tame` command line tool. 

## General usage

In general, the TAME commands takes a (or multiple) trajectory file as input.
For instance, the above command computes the self and distinct diffusion
coefficient, with the mean displacement correlation (MDC) of the atoms.

```
tame mdc traj.dump --seg 5000 --dt 5 --tags '1 1,1'
```

A handy feature is that all TAME recipes allows for the split of the
trajectories into segments of equal lengths, and perform the analysis on them
separately. For instance, the above command will output several output files:
`mcd_seg1.dat`, `mcd_seg2.dat`, ... the variance thereof serves as an error
estimation.

Most commands shared the options for input trajectories, while their specific
options can be found in the documentation pages, see below.

## List of commands

Below is a list of recipes implemented in TAME, with links to their respective
documentation pages.

| Command                      | Description                            |
|------------------------------|----------------------------------------|
| [`tame jacf`](jacf.md)       | current autocorrelation function       |
| [`tame mdc`](mdc.md)         | mean displacement correlation          |
| [`tame onsager`](onsager.md) | Onsager coefficients                   |
| [`tame pmsd`](pmsd.md)       | polarization mean squared displacement |
| [`tame persist`](persist.md) | persistent time                        |
| [`tame rdf`](rdf.md)         | radial distribution function           |
| [`tame sfac`](sfac.md)       | structure factor                       |

