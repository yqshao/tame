def cond_pmsd(dumps, pmsd_out, cond_out,
              unit, dt, max_dt, seg_dt, fit_min, fit_max, T,
              c_type, a_type, c_idx, a_idx):
    """Computing conductity using the Green-Kubo relation and the
    polarization mean squared displacement function

    Args:
        dumps: list lammps dump files.
        unit: dump file unit.
        dt: timestep of the trajectory [ps].
        max_dt: max correlation time to compute [ps].
        seg_dt: segment time to split the trajectory [ns].
        fit_min: 'minimal time for msd fit [ps]'.
        fit_max: 'maximal time for msd fit [ps]'.
        T: temperature in [K].
        c_type: cation element type (if None, must specify c_idx).
        a_type: anion element type (if None, must specify a_idx).
        c_idx: cation atom indices.
        a_idx: anion atom indices.
    """
    import numpy as np
    from tame.fit import linear_fit
    from tame.ops import tmsd, tavg, unwrap
    from tame.io.lammps import load_multi_dumps
    from tame.units import get_unit, ps, kB, e
    units = get_unit(unit)
    # Processing input
    n_segment = int(seg_dt*1e3/dt)
    n_cache = int(max_dt/dt)
    data = load_multi_dumps(dumps)
    if c_idx is None:
        assert c_type is not None, 'Must specify cation type or indices'
        c_idx = data['elems'] == c_type
    if a_idx is None:
        assert a_type is not None, 'Must specify anion type or indices'
        a_idx = data['elems'] == a_type
    # Operations
    P_c = unwrap(data['coord'][c_idx], data['cell'])
    P_a = unwrap(data['coord'][a_idx], data['cell'])
    P = np.sum(P_c - P_a, axis=0, keepdims=True)
    V = data['cell'][0]*data['cell'][1]*data['cell'][2]

    # Running
    count = 0
    result = []
    while True:
        count += 1
        pmsd = tmsd(P, n_cache)
        vavg = tavg(V)

        try:
            data.run(n_segment)
        except StopIteration:
            print('Reached end of trajectory')
            np.savetxt(f'{cond_out}.dat', result)
            print(f'Conductivity saved to {cond_out}.dat')
            break

        PMSD = pmsd.eval()
        VAVG = vavg.eval()
        TIME = np.arange(0, n_cache)*dt        
        slope, _ = linear_fit(TIME, PMSD, fit_min, fit_max)
        cond = 1/(6*kB*T*VAVG) * slope \
            /units.length**3 * e**2 * units.length**2 / ps
        result.append(cond)
        
        print(f'Segment {count} ({data.count+1} frames)'
              f': cond={cond:.2e}[S/m]', end='')
        if pmsd_out:
            outfile = f'{pmsd_out}_{count}.dat'
            np.savetxt(f'{outfile}', np.transpose([TIME, PMSD]))
            print(f', PMSD saved to "{outfile}.')
        else:
            print('.')


def set_parser(parser):
    from tame.recipes.utils import TypeDefaultFormat
    import argparse    
    parser.formatter_class=TypeDefaultFormat
    parser.description = 'Computing conductity using the '\
                         'current-current autocorrelation function.'
    parser.add_argument('dumps', metavar='dump', type=str, nargs='+',
                        help='dump files')
    parser.add_argument('--pmsd-out', type=str, default='pmsd',
                        help='polarization MSD data output')
    parser.add_argument('--cond-out', type=str, default='cond_pmsd',
                        help='conductivity data output')
    parser.add_argument('--unit', type=str, default='real',
                        help='dump file unit')
    parser.add_argument('--dt', type=float, default=0.5,
                        help='dump frequency in [ps]')
    parser.add_argument('--max-dt', type=float, default=20.0,
                        help='correlation time in [ps]')
    parser.add_argument('--seg-dt', type=float, default=5.0,
                        help='trajectory segment length in [ns]')
    parser.add_argument('--fit-min', type=float, default=5.0,
                        help='minimal time for msd fit [ps]')
    parser.add_argument('--fit-max', type=float, default=20.0,
                        help='maximal time for msd fit [ps]')    
    parser.add_argument('--temperature', type=float, default=273.15,
                        help='temperature in [K]')        
    parser.add_argument('--c-type', type=int, default=None,
                        help='cation atom type')
    parser.add_argument('--a-type', type=int, default=None,
                        help='anion atom type')
    parser.add_argument('--c-idx', type=int, default=None,
                        help='cation indices')
    parser.add_argument('--a-idx', type=int, default=None,
                        help='anion indices')       
    parser.set_defaults(func=lambda args: cond_pmsd(
        args.dumps, args.pmsd_out, args.cond_out,
        args.unit, args.dt, args.max_dt, args.seg_dt,
        args.fit_min, args.fit_max, args.temperature,
        args.c_type, args.a_type, args.c_idx, args.a_idx))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
