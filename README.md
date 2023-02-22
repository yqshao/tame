# TAME: Trajectory Analysis Made Easy

## Introduction

This is a **work in progress** set of tools to post process MD trajectories,
mainly aiming for calculating transport properties of electrolyte systems.

## Installation

**with pip**

``` shell
   pip install git+https://github.com/Teoroo-CMC/tame.git
```

**with singularity**

``` shell
   singularity build (--remote) tame.sif docker:yqshao/tame:latest
```

## Using a recipe

Recipes are organzied in subcommands of the mdppp command, for example, the
following command computes the Green-Kubo conductivity taking element 1 and 2 as
cation and anion species:

``` shell
tame diff mcd traj.dump --seg 5000 --dt 5 --tag '1 1,1'
```
