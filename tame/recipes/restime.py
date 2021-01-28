#!/usr/bin/env python3

def restime(dumps, top, tcf_out, unit, dt, max_dt, seg_dt, tags):
    """Computes the time correlation function

    Args:
        dumps: list lammps dump files.
        unit: dump file unit.
        tcf_out: naming for residence time correlation function.
        dt: timestep of the trajectory [ps].
        max_dt: max correlation time to compute [ps].
        seg_dt: segment time to split the trajectory [ns].
        tags: list of residence times to compute
    """
    import numpy as np
    from tame.ops import tcache, tavg, tsurvive, unwrap
    from tame.io import load_traj
    from tame.units import get_unit, ps, kB, e
    units = get_unit(unit)
    # Processing input
    n_segment = int(seg_dt*1e3/dt)
    n_cache = int(max_dt/dt)
    data = load_traj(dumps, top)

    coord = unwrap(data['coord'], data['cell'])
    cell = data['cell'][None, None,:]

    # Defining the correlations to evaluate
    def make_corrs():
        corrs = {}
        all_dr = {}
        dists = {}
        in_cache = {}
        for tag in tags:
            t1, t2 = map(int, tag.split(':')[0].split(','))
            method = tag.split(':')[1]
            rc = [*map(float, tag.split(':')[2].split(','))]
            # All the pairs we are going to care
            idx_i = np.where((data['elems']==t1).eval())[0]
            idx_j = np.where((data['elems']==t2).eval())[0]
            r_i = coord[idx_i]; r_j = coord[idx_j]
            if (t1, t2) not in dists:
                r_ij = r_i[:,None,:] - r_j[None, :, :]
                r_ij = r_ij - np.rint(r_ij/cell)*cell
                dist = np.sqrt(np.sum(r_ij**2, axis=2))
                dists[(t1, t2)] = dist
            birth = dists[(t1, t2)]<=rc[0]
            death = dists[(t1, t2)]>rc[int(method=='SSP')]
            if t1==t2:
                birth = birth * (1-np.eye(len(idx_i))[None, :, :])
            if method == 'IMM':
                n_death = int(rc[1]/dt)
                death = np.logical_and.accumulate(
                    tcache(death, n_death), axis=0)
                corrs[tag] = tavg(np.logical_and.accumulate(
                    tcache(tsurvive(birth, death), n_cache), axis=0)
                                  .sum(axis=2).mean(axis=1))
            elif method == 'SSP':
                corrs[tag] = tavg(np.logical_and(
                    tcache(birth, n_cache), np.logical_and.accumulate(
                        tcache(~death, n_cache), axis=0))
                                  .sum(axis=2).mean(axis=1))
        return corrs

    # Running
    count = 0
    result = []
    while True:
        corrs = make_corrs()
        keys = list(corrs.keys())
        labels = [f'D_{{{k}}}' for k in list(corrs.keys())]
        count += 1

        try:
            data.run(n_segment)
        except StopIteration:
            print('Reached end of trajectory')
            break

        TIME = np.arange(0, n_cache)*dt
        CORRS = [corrs[k].eval() for k in keys]
        outfile = f'{tcf_out}_{count}.dat'
        np.savetxt(f'{outfile}', np.transpose([TIME]+CORRS))
        print(f'Correlation functions saved to "{outfile}.')

def set_parser(parser):
    from tame.recipes.utils import TypeDefaultFormat
    import argparse
    parser.formatter_class = TypeDefaultFormat
    parser.description = """
Three types of correation functions supported corresponding to three different
definitions of residence time. More details can be found in: Lagge and
Hynes, J. Phys. Chem. B, 2008, 112, 7697-7701.

RF (reactive flux): e.g. `3,4:SF:3.0`
-------------------------------------
The reactive flux method, note that here only we consider only the first passage
time, thus RF is equivatlent to IMM with t*=0 or SSP with a single cutoff.

IMM (Impey, Madden, McDonald): e.g. `3,4:IMM:3.0,1`
---------------------------------------------------
The pairs are defined as two ions that doees not depart for more than certain
time takes a cutoff distance and a time (in ps) as argument.

SSP (stable state picture): e.g. `3,4:SSP:2.8,3.2`
--------------------------------------------------
The correlation function is given by the probability of originating from a
stable reactant and not reaching a stable product, takes two cutoff distances as
arguments.
"""

    parser.add_argument('dumps', metavar='dump', type=str, nargs='+',
                        help='dump files')
    parser.add_argument('--top', type=str, default='',
                        help='topology file')
    parser.add_argument('--tcf-out', type=str, default='restcf',
                        help='residence time correlation fundtion output')
    parser.add_argument('--unit', type=str, default='real',
                        help='dump file unit')
    parser.add_argument('--dt', type=float, default=0.5,
                        help='dump frequency in [ps]')
    parser.add_argument('--max-dt', type=float, default=20.0,
                        help='correlation time in [ps]')
    parser.add_argument('--seg-dt', type=float, default=5.0,
                        help='dump frequency in [ps]')
    parser.add_argument('--tags', type=str, nargs='+', default=[],
                        help='types of residence times')

    parser.set_defaults(func=lambda args: restime(
        args.dumps, args.top, args.tcf_out, args.unit,
        args.dt, args.max_dt, args.seg_dt, args.tags))

def main():
    import argparse
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
