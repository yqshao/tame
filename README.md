# TAME: Trajectory Analysis Made Easy

## Introduction

Tame is a set of tools to post process MD trajectories, mainly aiming for
calculating transport properties of electrolyte systems from time correlation
functions. File IO is through the [chemfiles] library.

[chemfiles]: https://chemfiles.org/

## Installation

**with pip**

``` shell
   pip install git+https://github.com/yqshao/tame.git
```

**with singularity**

``` shell
   singularity build (--remote) tame.sif docker:yqshao/tame:latest
```

## Using a recipe

Recipes are organized in subcommands of the `tame` command, for example, the
following command computes the Green-Kubo conductivity taking element 1 and 2 as
cation and anion species:

``` shell
tame diff mcd traj.dump --seg 5000 --dt 5 --tag '1 1,1'
```
