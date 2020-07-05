# Basic concepts

## Trajectory and FrameArray

### Trajectory

In TAME, a trajectory is a continuous stream of atoms, the `Trajectory` object
can be loaded from one or a ordered list of trajectory files from MD programs, 
e.g:
    
```Python
from mdppp.io import load_traj
traj = load_traj('prod.dump')
```

### FrameArray


A `FrameArray` means a numpy-array-like data structure which is updated whenever
a `Trajectory` (called its parent) is updated. The ``FrameArray`` behaves much
like a numpy array, with the exception that its mathematical operations produces
`FrameArray`s (also updated according to the udpated value of itself every
frame).


```Python
r0 = traj.coord(0)
sin_r0 = np.sin(r0)
for i in range(3):
    data.run(1)
    print(r0, sin_r0)
```

### Under the hood

When we calculate from s `FrameArray`, a computation graph is implicitly
constructed in the `Trajectroy` object. Some for example, the `coord` method
retrieves the coordinate of some atoms from the trajectory. Whenever
`traj.run()` is called, all `FrameArray`s are updated according to the order in
which they are created.

##  Time-dependent functions

### Using a time-depend function

TAME computes structure properties one-frame at a time. However, we are often
more interested in some property that depend on previous frames, e.g. the time
correlation functions.

The most common time-dependent function we'd like to use is the `tavg` function,
which gives the cumulative average value of a `FrameArray`. Here, we show how
one get the ensemble averaged radial distribution function (RDF).

```Python
r_grid, rdf = tavg(traj.rdf('O', 'O', 2, 5, 100))
traj.run(100)
plt.plot(r_grid, rdf)
traj.run(900)
plt.plot(r_grid, rdf)
```

One see how the RDF converges with longer trajectory.

### Available time functions

All implemented time-dependent functions are in the `tame.time` module. See the
API documentation for more detail.

### Writing a custom time function

Unlike the automatically computed static `FrameArray`, each time-dependent
function requires different operations when a new frame is loaded. Therefore,
all custom time function was implemented as a function producing a special
subclass of `FrameArray`. For instance, a simplified version of the `tavg`
is shown below:

```Python
class TAVG(FrameArray):
    def __init__(self, var):
        var.parent.derived.append(self)
        self.parent = var.parent
        self.var = var
        self.cumsum = 0
        self.count = 0
        self.update()

    def update(self):
        new_val = self.var.numpy()
        self.count += 1

    def numpy(self):
        return self.cumsum/self.count

def tavg(var, dropnan='partial'):
    return TAVG(var)
```

A `TAVG` object keeps track of the cumulative sum of the variable and its mean
value when a new frame is loaded, only when we call the `.numpy()` method is the
actually mean value calculated.
