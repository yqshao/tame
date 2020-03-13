"""This module implements loaders for lammps trajectry and log files."""

import numpy as np

def _dump_generator(fname, resort=False):
    f = open(fname)
    [next(f) for _ in range(8)]
    l = f.readline()
    fmt = l.split()[2:]
    idx_coord = [fmt.index(key) for key in ['x','y','z']]
    idx_speed = [fmt.index(key) for key in ['vx','vy','vz']]
    idx_elem = fmt.index('type')
    f.close()
    f = open(fname)
    [next(f) for _ in range(3)]
    natoms = int(f.readline())
    def get_cell(line): return float(line.split()[1]) - float(line.split()[0])
    count = 0
    while True:
        line = f.readline()
        if line == '':
            break
        if line.startswith('ITEM: BOX'):
            cell = np.array([get_cell(f.readline()) for i in range(3)])
            next(f)
            count += 1
            coord, speed, elems = [], [], []
            for i in range(natoms):
                line = f.readline().split()
                coord.append([line[idx] for idx in idx_coord])
                speed.append([line[idx] for idx in idx_speed])
                elems.append(line[idx_elem])
            elems = np.array(elems, np.int32)
            data = {
                'coord': np.array(coord, np.float64),
                'speed': np.array(speed, np.float64),
                'elems': elems,
                'cell': cell}
            yield data


def _multi_dump_generator(flist, *args, **kwargs):
    """Generator for multiple dump files"""
    for i, dump in enumerate(flist):
        for j, data in enumerate(_dump_generator(dump)):
            if i == 0 or j > 0:
                yield data


def load_multi_dumps(flist, resort=False):
    """Read a list of file, gives a FrameData

    Expected lammps dump format: "ATOMS id type x y z vx vy vz"
    Atoms should be sorted with index

    Args: 
        flist: name of dump file

    Returns:
        A FrameData containing:
        - coord: coordinate array, size=(natoms, 3)
        - speed: velocity array, size=(natoms, 3)
        - elems: element array, size=(natoms)
        - cell: cell array, size=(3, 3)
    """
    from mdppp import FrameData
    return FrameData(_multi_dump_generator(flist))
