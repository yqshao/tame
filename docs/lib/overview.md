# TAME Library

TAME can be used as a python library to compose your analysis workflows. This
part of the documentation introduces different components of TAME, and how you
can extend TAME to add your analysis.

## Concepts

**Analysis** Most of TAME's recipes may be represented as an `Analysis` object,
which completely defines a sequence of "operations" including computation,
caching and outputting. "Running" an analysis means going through a trajectory,
and updating the `TArrays` by executing each `Operation` in sequence at every
step.

**TArrays** In the analyses, data are represented with so-called `TArray`s in
TAME. TArrays look and behave like numpy arrays; you can run perform simple
arithmetic operations and numpy
[ufunc](https://numpy.org/doc/stable/reference/ufuncs.html) on them, which in
turn produces TArrays instead of numpy arrays. Unlike numpy arrays, TArrays
represent the state of analyses and are therefore updated whenever the analyses
proceed.

**Operations** The mathematical operations are simple cases of the so-called
`Operations` in TAME. Besides them, some operations accumulate results within an
analysis, produce the visualization (for GUI), or write out the results (for CLI
recipes).

## Modules 

The TAME library is organized as follows:

- `tame.core`: Core implementations of the `TArray`, `Analysis`;
- `tame.ops`: Specific `Operations` for molecular analyses;
- `tame.io`: Loading of trajectories (based on MDAnalysis).