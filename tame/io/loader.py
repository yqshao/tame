"""`tame.io.loader`

The loader module implments the functions to create the `Trajectory` object
from trajectory files, which includes:

- `load_traj`: A univeral loader for trajectries from different MD codes, see
  the function documentation for supported formats.
- `load_random_walk`: A loader generating a Trajectory of random walking
  particles, this loader is mainly for the purpose of test.

"""
def guess_format(fname):
    import re
    if re.search(r'\.dump(\.xz|\.gz)?$', fname):
        return 'lammps-dump'
    if fname.endswith('.xtc'):
        return 'gromacs-trr'

def load_traj(trajectory, topology=None, format='auto'):
    """ Trajectory file loader

    By default, the file format is guessed from the file name of the first
    trajectory file, this behavior can be tweaked with the format option.

    Supported trajectory formats

    | Format         | Trajectory Suffix | Description            |
    |----------------|-------------------|------------------------|
    | `lammps-dump`  | `.dump`           | LAMMPS trajectory file |
    | `gromacs-trr`  | `.trr`            | GROMACS trajectory     |

    Args:
        trajectory (str or list): (list of) trajectory files
        topology (str): topology file
        format (str): file format

    Returns:
        (Trajectory): TAME Trajectory object
    """
    from tame import io

    if isinstance(trajectory, str):
        trajectory = [trajectory]

    if format == 'auto':
        format = guess_format(trajectory[0])

    loader = {
        'lammps-dump': lambda: io.load_lammps_dump(trajectory),
        'gromacs-trr': lambda: io.load_gromacs_trr(trajectory, topology)
    }[format]
    traj = loader()

    return traj


def load_random_walk(n_particles, pbc=False, cell=None):
    """Generating a trajectory from random walk

    TO BE IMPLEMENTED.
    """
    pass
