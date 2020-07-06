def cond_jacf(dumps, jacf_out, cond_out,
              unit, dt, max_dt, seg_dt, T,
              c_type, a_type, c_idx, a_idx):
    """Computing conductity using the Green-Kubo relation and the
    current-current autocorrelation function

    Unit conversion:
        JACF integral is given in speed*2 * ps.
        V is given in lammps volume unit.
        T is given in K.
        Returning conductivity in SI [S/m].

    Args:
        dumps: list lammps dump files.
        unit: dump file unit.
        dt: timestep of the trajectory [ps].
        max_dt: max correlation time to compute [ps].
        seg_dt: segment time to split the trajectory [ns].
        T: temperature in [K].
        c_type: cation element type (if None, must specify c_idx).
        a_type: anion element type (if None, must specify a_idx).
        c_idx: cation atom indices.
        a_idx: anion atom indices.
    """
    import numpy as np
    from tame.ops import tacf, tavg
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
    J_c = data['speed'][c_idx]
    J_a = data['speed'][a_idx]
    J = np.sum(J_c - J_a, axis=0, keepdims=True)
    V = data['cell'][0]*data['cell'][1]*data['cell'][2]

    # Running
    count = 0
    result = []
    while True:
        count += 1
        jacf = tacf(J, n_cache)
        vavg = tavg(V)

        try:
            data.run(n_segment)
        except StopIteration:
            print('Reached end of trajectory')
            np.savetxt(f'{cond_out}.dat', result)
            print(f'Conductivity saved to {cond_out}.dat')
            break

        JACF = jacf.eval()
        VAVG = vavg.eval()
        cond = 1/(3*kB*T*VAVG)*np.trapz(JACF, dx=dt)\
            /units.length**3 * e**2 * units.velocity**2 * ps
        result.append(cond)

        print(f'Segment {count} ({data.count+1} frames)'
              f': cond={cond:.2e}[S/m]', end='')
        if jacf_out:
            TIME = np.arange(0, n_cache)*dt
            outfile = f'{jacf_out}_{count}.dat'
            np.savetxt(f'{outfile}', np.transpose([TIME, JACF]))
            print(f', JACF saved to "{outfile}.')
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
    parser.add_argument('--jacf-out', type=str, default='jacf',
                        help='current autocorrelation data output')
    parser.add_argument('--cond-out', type=str, default='cond_jacf',
                        help='conductivity data output')
    parser.add_argument('--unit', type=str, default='real',
                        help='dump file unit')
    parser.add_argument('--dt', type=float, default=0.5,
                        help='dump frequency in [ps]')
    parser.add_argument('--max-dt', type=float, default=20.0,
                        help='correlation time in [ps]')
    parser.add_argument('--seg-dt', type=float, default=5.0,
                        help='trajectory segment length in [ns]')
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
    parser.set_defaults(func=lambda args: cond_jacf(
        args.dumps, args.jacf_out, args.cond_out,
        args.unit, args.dt, args.max_dt, args.seg_dt, args.temperature,
        args.c_type, args.a_type, args.c_idx, args.a_idx))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
