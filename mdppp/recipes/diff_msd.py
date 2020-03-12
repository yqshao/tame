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
    from mdppp.fit import linear_fit
    from mdppp.ops import tmsd, tavg, unwrap
    from mdppp.io.lammps import load_multi_dumps
    from mdppp.units import get_unit, ps, kB, e
    units = get_unit(unit)
    # Processing input
    n_segment = int(seg_dt*1e3/dt)
    n_cache = int(max_dt/dt)
    data = load_multi_dumps(dumps)

    coord = unwrap(data['coord'], data['cell'])
    # Defining the correlations to evaluate
    for t in types:
        if '-' not in t: # just simple self diffusion
            corrs[t] = tmsd(coord[data['elems']==int(t)])
        else:
            t1, t2 = map(int, t.split('-'))
            if t1==t2: # same species
                r_i = coord[data['elems']==t1]
                dr_i = tcache(r_i, n_cache) - r_i[None, :]
                dr_j = np.sum(dr_i, axis=0, keepdims=True)-dr_i
                corr = dr_i[None,:,:]*tcache(dr_j, n_cache)
                corr = np.mean(np.sum(corr, axis=1), axis=2)
            else: # Different species
                dr_i = np.sum(coord[data['elems']==t1], axis=0, keepdims=True)
                dr_j = np.sum(coord[data['elems']==t2], axis=0, keepdims=True)
                corr = dr_i[None,:,:]*tcache(dr_j, n_cache)
                corr = np.mean(np.sum(corr, axis=1), axis=2)

        corrs[t] = tavg[corr]

    # Running
    count = 0
    result = []
    while True:
        count += 1
        corrs = tmsd(P, n_cache)

        try:
            data.run(n_segment)
        except StopIteration:
            print('Reached end of trajectory')
            np.savetxt(f'{diff_out}.dat', result)
            print(f'Diffusion coefficients saved to {diff_out}.dat')
            break

        CORRS = [corr.eval() k in keys()]
        diffs = {1/6*linear_fit(TIME, CORR, fit_min, fit_max)
                 for CORR in CORRS}
        result.append(diffs)

        print(f'Segment {count} ({data.count+1} frames)'
              ' '.join([f':, {key}={diff:.2e}'
                        for key, diff in zip(keys,diffs)], end='[Ang^2/ps]')
        if pmsd_out:
            outfile = f'{corr_out}_{count}.dat'
            np.savetxt(f'{outfile}', np.transpose([TIME]+CORRS))
            print(f', Correlation functions saved to "{outfile}.')
        else:
            print('.')

def set_parser(parser):
    from mdppp.recipes.utils import TypeDefaultFormat
    import argparse
    parser.formatter_class=TypeDefaultFormat
    parser.description = 'Computing various diffusion coefficients'

    parser.add_argument('dumps', metavar='dump', type=str, nargs='+',
                        help='dump files')
    parser.add_argument('--corr-out', type=str, default='pmsd',
                        help='correlation function output')
    parser.add_argument('--coeff-out', type=str, default='cond_pmsd',
                        help='coefficient output')
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
    parser.add_argument('--types', type=int, nargs='+',
                        help='types of particles')
    parser.add_argument('--by-bond', type=str, default=0,
                        help='divide the distinct diffusion by bond/nonbond')
    parser.add_argument('--by-cutoff', type=str, default=0,
                        help='divide the distinct diffusion by in/out')

    parser.set_defaults(func=lambda args: diff_msd(
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
