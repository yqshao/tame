"""`tame.io.loader`

The loader module implments the functions to create the `Trajectory` object
from trajectory files, which includes:

- `load_traj`: A univeral loader for trajectries from different MD codes, see
  the function documentation for supported formats.
- `load_random_walk`: A loader generating a Trajectory of random walking
  particles, this loader is mainly for the purpose of test.

"""

def load_traj(trajectory, topology=None, timestep=1, format=None):
    """ Trajectory file loader

    By default, the file format is guessed from the file name of the first
    trajectory file, this behavior can be tweaked with the format option.

    Supported trajectory formats

    | Format    | Trajectory Suffix | Description            |
    |-----------|-------------------|------------------------|
    | `lammps`  | `.dump`           | LAMMPS trajectory file |
    | `gromacs` | `.trr`            | GROMACS trajectory     |

    Args:
        trajectory (str or list): (list of) trajectory files
        topology (str): topology file
        timestep (float): timestep of the trajectory
        format (str): file format

    Returns:
        (Trajectory): TAME Trajectory object
    """
    return traj


def load_random_walk(n_particles, pbc=False, cell=None):
    """Generating a trajectory from random walk

    TO BE IMPLEMENTED.
    """
    pass
