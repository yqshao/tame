"""This recipe computes the ionic conductivity using the Green-Kubo relation."""

def cond_jacf(dumps, jacf_out, cond_out,
              unit, dt, max_dt, seg_dt,
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
        c_type: cation element type (if None, must specify c_idx).
        a_type: anion element type (if None, must specify a_idx).
        c_idx: cation atom indices.
        a_idx: anion atom indices.
    """
    import numpy as np
    from mdppp.op import tacf
    from mdppp.io.lammps import load_multi_dumps
    # from mdppp.units import get_unit, ps, e    
    # Processing input
    n_segment = seg_dt//dt
    n_cache = max_dt//dt    
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
    V = data['cell']**3
    # Running
    count = 0
    result = []
    while True:
        count += 1
        jacf = tacf(J, n_cache)
        vavg = tavg(V)

        if not data.run(n_segment):
            print('Reached end of trajectory')
            np.savetxt(f'{cond_out}.dat', result)
            print(f'G-K conductivity saved to {cond_out}.dat')
            break

        JACF = jacf.eval()
        V = vavg.eval()
        cond = np.trapz(JACF, dx=dt) #* unit.speed**2*ps/unit.volume
        result.append(cond)

        print(f'Segment {count} ({data.count} frames)'
              f': cond={cond:.2e}[S/m]', end='')
        if jacf_log:
            TIME = np.arange(0, cache_size)*dt
            outfile = f'{jacf_out}_{count}.dat'
            np.savetxt(f'{outfile}', [TIME, JACF])
            print(f', JACF saved to "{outfile}."')
        else:
            print('.')


def set_parser(parser):
    from mdppp.recipe.utils import TypeDefaultFormat
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
        args.unit, args.dt, args.max_dt, args.seg_dt,
        args.c_type, args.a_type, args.c_idx, args.a_idx))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
