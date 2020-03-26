def diff_dcf(dumps, diff_out, corr_out,
              unit, dt, max_dt, seg_dt, fit_min, fit_max, T,
              tags):
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
    from mdppp.ops import tmsd, tcache, tavg, unwrap
    from mdppp.io.lammps import load_multi_dumps
    from mdppp.units import get_unit, ps, kB, e
    units = get_unit(unit)
    # Processing input
    n_segment = int(seg_dt*1e3/dt)
    n_cache = int(max_dt/dt)
    data = load_multi_dumps(dumps)

    coord = unwrap(data['coord'], data['cell'])

    # Defining the correlations to evaluate
    corrs = {}
    all_dr = {}
    def get_dr(t):
        if t not in all_dr.keys():
            r_t = coord[data['elems'].eval()==t]
            all_dr[t] = tcache(r_t, n_cache) - r_t[None, :, :]
        return all_dr[t]

    dists = {}
    cell = data['cell'][None, None,:]
    for tag in tags:
        if ',' not in tag: # just simple self diffusion
            corr = np.mean(np.sum(get_dr(int(tag))**2, axis=2), axis=1)
            corrs[tag] = tavg(corr)
        elif ':' not in tag: #distinct diffusion
            t1, t2 = map(int, tag.split(','))
            if t1==t2: # same species
                dr_i = get_dr(t1)
                dr_j = np.sum(dr_i, axis=1, keepdims=True) - dr_i
                corr = np.mean(np.sum(dr_i*dr_j, axis=2), axis=1)
            else: # Different species
                dr_i = np.mean(get_dr(t1), axis=1)
                dr_j = np.sum(get_dr(t2), axis=1)
                corr = np.sum(dr_i*dr_j, axis=1)
            corrs[tag] = tavg(corr)
        elif ':C' in tag: # Handel the cutoffs here
            ts, r_c = tag.split(':C')
            r_c = float(r_c)
            t1, t2 = map(int, ts.split(','))
            # All the pairs we are going to care
            idx_i = np.where((data['elems']==t1).eval())[0]
            idx_j = np.where((data['elems']==t2).eval())[0]
            r_i = coord[idx_i]; r_j = coord[idx_j]
            if (t1, t2) not in dists:
                r_ij = r_i[:,None,:] - r_j[None, :, :]
                r_ij = r_ij - np.rint(r_ij/cell)*cell
                dist = np.sqrt(np.sum(r_ij**2, axis=2))
                dists[(t1, t2)] = dist
            # "Memory" of in/out
            in_cache = tcache(dists[(t1, t2)]<=r_c, n_cache)
            # Criterion for in/out
            cut_in = np.logical_and.accumulate(in_cache, 0)
            cut_out =  np.logical_not(cut_in)
            if t1==t2:
                cut_in = cut_in * (1-np.eye(len(idx_i)))
            # Add the corresponding j_in/out for each i first
            dr_i = get_dr(t1); dr_j = get_dr(t2)
            dr_j_in  = np.sum(dr_j[:,None,:,:]* cut_in[:,:,:,None],axis=2)
            dr_j_out = np.sum(dr_j[:,None,:,:]*cut_out[:,:,:,None],axis=2)
            corrs[tag+f'_in'] = tavg(np.mean(np.sum(dr_i*dr_j_in, axis=2), axis=1))
            corrs[tag+f'_out'] = tavg(np.mean(np.sum(dr_i*dr_j_out, axis=2), axis=1))
        else:
            raise NotImplementedError(f"We can not handel the tag {tag}, sorry dude.")

    keys = list(corrs.keys())
    labels = [f'D_{{{k}}}' for k in list(corrs.keys())]

    # Running
    count = 0
    result = []
    while True:
        count += 1

        try:
            data.run(n_segment)
        except StopIteration:
            print('Reached end of trajectory')
            np.savetxt(f'{diff_out}.dat', result, header=' '.join(labels))
            print(f'Diffusion coefficients saved to {diff_out}.dat')
            break

        TIME = np.arange(0, n_cache)*dt
        CORRS = [corrs[k].eval() for k in keys]
        diffs = [1/6*linear_fit(TIME, CORR, fit_min, fit_max)[0]
                 for CORR in CORRS]
        result.append(diffs)

        print(f'Segment {count} ({data.count+1} frames)'+
              ''.join([f', {label}={diff:.2e}'
                        for label, diff in zip(labels,diffs)]), end=' [Ang^2/ps]')
        if corr_out:
            outfile = f'{corr_out}_{count}.dat'
            np.savetxt(f'{outfile}', np.transpose([TIME]+CORRS))
            print(f', Correlation functions saved to "{outfile}.')
        else:
            print('.')

def set_parser(parser):
    from mdppp.recipes.utils import TypeDefaultFormat
    import argparse
    parser.formatter_class = TypeDefaultFormat
    parser.description = """Computing various diffusion coefficients
  self diffusion coefficients (D^s) can be specified using their type
  disctinct diffusion coefficients (D^d) is specified using e.g. 3,3

  you can further split the D^d by a cutoff radius using e.g. 3,3:C3.5"""

    parser.add_argument('dumps', metavar='dump', type=str, nargs='+',
                        help='dump files')
    parser.add_argument('--diff-out', type=str, default='diff',
                        help='dffusion coefficients output')
    parser.add_argument('--corr-out', type=str, default='corr',
                        help='correlation functions output')
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
    parser.add_argument('--tags', type=str, nargs='+', default=[],
                        help='types of diffusion coefficient to calculate')

    parser.set_defaults(func=lambda args: diff_dcf(
        args.dumps, args.diff_out, args.corr_out,
        args.unit, args.dt, args.max_dt, args.seg_dt,
        args.fit_min, args.fit_max, args.temperature,
        args.tags))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
