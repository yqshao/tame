"""`tame.io.writer`

The writer module implements `write_traj()` function with generates a trajectory
file from `FrameArray`s. See the function documenatation for supported formats.
"""

def write_traj(fname, coord, freq=1, types=None, format=None):
    """Trajectory writer

    TO BE IMPLEMENTED

    The `write_traj` function generates a `FrameArray` subclass whose `update`
    method writes to a trajectory file.

    Supported trajectory formats:
    - None

    Args:
        fname (str): file name of the output trajectroy
        freq (int): dumping frequency in frames
        types (list): list of atom types

    Returns:
        (FrameWriter): subclass of FrameArray which dumps the trajectory
    """
    pass
