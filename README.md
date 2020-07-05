# Trajectory Analysis Made Easy

## Introduction

This is a **work in progress** set of tools to post process MD trajectories,
mainly aiming for calculating transport properties of electrolyte systems.

## Installation

``` shell
   pip install git+https://github.com/yqshao/mdppp.git
```

## Using a recipe

Recipes are organzied in subcommands of the mdppp command, for example, the
following command computes the Green-Kubo conductivity taking element 1 and 2 as
cation and anion species:


``` shell
   mdppp cond jacf prod.dump --c-type 1 --a-type 2
```

