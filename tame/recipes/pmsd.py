#!/usr/bin/env python3
import click
import tame
from tame.recipes.utils import load_traj_seg

@click.command(name='pmsd',
               options_metavar='[options]',
               short_help='polarization mean squared displacement')
@load_traj_seg # general input handler
@click.option('--max-dt',       metavar='', default=30.0,   show_default=True)
@click.option('-c', '--c-type', metavar='', default=3,    show_default=True)
@click.option('-a', '--a-type', metavar='', default=4,    show_default=True)
@click.option('-T', '--temp',   metavar='', default=297.3,  show_default=True)
@click.option('--fit-min',      metavar='', default=5.0,    show_default=True)
@click.option('--fit-max',      metavar='', default=20.0,   show_default=True)
@click.option('--pmsd-out',     metavar='', default='pmsd', show_default=True)
@click.option('--cond-out',     metavar='', default='cond', show_default=True)
def pmsd_cmd(seg, dt, # provided by load_traj_seg
             temp, c_type, a_type, max_dt, fit_min, fit_max, pmsd_out, cond_out):
    """Computing conductity using the Green-Kubo and the Einstein-relation, from the
    polarization mean squared displacement (PMSD).

    See the documentation for detailed descriptions of the command:
    https://Teoroo-CMC.github.io/tame/latest/recipe/pmsd

    """
    import numpy as np
    from tame.ops import tmsd, tavg, unwrap
    from tame.fit import linear_fit
    from tame.constants import units, kB, e

    n_cache = int(max_dt/dt)
    c_idx = seg['elems'] == c_type
    a_idx = seg['elems'] == a_type
    P_c = unwrap(seg['coord'][c_idx], seg['cell'])
    P_a = unwrap(seg['coord'][a_idx], seg['cell'])
    P = np.sum(P_c - P_a, axis=0, keepdims=True)
    pmsd = tmsd(P, n_cache)
    V = seg['cell'][0]*seg['cell'][1]*seg['cell'][2]
    Vavg = tavg(V)

    seg.run()
    PMSD = pmsd.eval()
    VAVG = Vavg.eval()
    TIME = np.arange(0, n_cache)*dt
    slope, _ = linear_fit(TIME, PMSD, fit_min, fit_max)
    cond = 1/(6*kB*temp*VAVG) * slope \
        /units.length**3 * e**2 * units.length**2 / units.time
    return {pmsd_out: {'t':TIME, 'PMSD': PMSD}, cond_out: cond}
