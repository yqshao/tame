# TAME Libarary

This section describes the basic concepts to use TAME as trajectory analysis
library to implement your own analysis. 

## Example analysis

Below is a short demo of how an analysis is composed and executed in TAME. The
code block that defines hydroxyl ion as the oxygens boned to only one hydrogen,
and survives at least 1 ps. The analysis then computes the radial distribution
functions for the hydroxyl ions.

```Python
from tame.io import load_demo
import tame.ops as ops
import tame.units as u

H = ops.sel_elem('H')
O = ops.sel_elem('O')
CN_O = ops.assign(H, O)
OH = ops.by_bool(CN_O==1)
OH = ops.t_persist(OH, 1*u.ps)
analysis = ops.collect(
    ops.rdf(OH, H, fout='rdf_OH'),
    ops.rdf(OH, OH, fout='rdf_OO')
)

traj = load_demo('water.xyz')
analysis.run(traj, to=100*u.ps)
```

## How this works

The core data structures of TAME are `Trajectory` and `Operation`. The
`Trajectory`s hold the data and metadata of a MD trajectory, while an
`Operation` defines a function that acts like an function of `Trajectory`.
Unlike a pure function, an `Operation` is stateful, which means it can be
history-dependent, such as a time average or a time correlation. As shown above,
the `Operation` acts like a numpy array upon arithmetic operaitons, and yields
new `Operation` objects.

## Further reading

The `Trajectory` and `Operation` provides a flexible way to define your MD
analysis with complex and structural dependence. See the **Examples** section
for some runnable Jupyter Notebooks. The **Data types** and **Module** sections
contain more detailed API documentation.

