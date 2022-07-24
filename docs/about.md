# About TAME

TAME (Trajectory Analysis Made Easy) was built by Yunqi Shao at Uppsala
Unversity.

The program was inspired by Matti Hellström's (in-house) analysis code for MD
trajectories, which has a much flexible way to define groups of atoms, and made
available a wide range of structural and dynamic properties (lifetime of
hydrogen bonds, diffusion of reactive species like hydroxyl ions, to name a
few). Despite its high performance and flexibility, the syntax was a bit tricky
to grasp and the program was written in C++, both made it hard to extend the
code.

Such difficulties hindered us when we start develop new gauges to the MD
trajectories, for instance, the cutoff-separated distinct diffusion
coefficients. We therefore developed a new set of tools such that the analysis
can be written in the more accessible language of Numpy arrays. The result was
TAME. The first prototype of the code was developed by Yunqi Shao, with the help
from Supho Phunnarungsi and Harish Gudla.

## Similar packages

Below is a list of packages that provides similar functionalities as TAME.

- [LiquidLib](https://z-laboratory.github.io/LiquidLib/) by Z Lab 
- [MDAnalysis](https://www.mdanalysis.org/)
- [MDTraj](https://www.mdtraj.org)
- [TRAVIS](http://www.travis-analyzer.de/)
