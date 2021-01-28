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
    from tame.ops import tavg, tpairsurvive
    from tame.io import load_traj
    from tame.units import get_unit, ps, kB, e
    units = get_unit(unit)
    # Processing input
    n_segment = int(seg_dt*1e3/dt)
    n_cache = int(max_dt/dt)
    data = load_traj(dumps, top)
    d_cache = {}
    make_corrs = lambda: {tag: tavg(np.mean(np.sum(tpairsurvive(
        data, tag, n_cache, d_cache), axis=2), axis=1)) for tag in tags}

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
    parser.description = """ See tame.time.tpairsurvive for details about the definition of survival
    probability."""

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
