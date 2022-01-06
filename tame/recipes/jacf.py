#!/usr/bin/env python3
import click
import tame
from tame.recipes.utils import load_traj_seg

@click.command(name='jacf',
               options_metavar='[options]',
               short_help='current autocorrelation function')
@load_traj_seg # general input handler
@click.option('--max-dt',       metavar='', default=30.0,   show_default=True)
@click.option('-c', '--c-type', metavar='', default='3',    show_default=True)
@click.option('-a', '--a-type', metavar='', default='4',    show_default=True)
@click.option('-T', '--temp',   metavar='', default=297.3,  show_default=True)
@click.option('--jacf-out',     metavar='', default='jacf', show_default=True)
@click.option('--cond-out',     metavar='', default='cond', show_default=True)
def jacf_cmd(seg, dt, # provided by load_traj
             c_type, a_type, max_dt,  T, jacf_out, cond_out):
    """Computing conductity using the Green-Kubo relation and the current-current
    autocorrelation function.

    See the documentation for detailed descriptions of the command:
    https://Teoroo-CMC.github.io/tame/latest/recipe/jacf
    """
    import numpy as np
    from numpy import trapz
    from tame.ops import tacf, tavg
    from tame.constants import units, kB, e

    n_cache = int(max_dt/dt)
    c_idx = seg['elems'] == c_type
    a_idx = seg['elems'] == a_type
    J_c = seg['speed'][c_idx]
    J_a = seg['speed'][a_idx]
    J = np.sum(J_c - J_a, axis=0, keepdims=True)
    V = seg['cell'][0]*seg['cell'][1]*seg['cell'][2]
    Vavg = tavg(V)

    seg.run()
    JACF = jacf.eval()
    VAVG = Vavg.eval()
    TIME = np.arange(0, n_cache)*dt
    COND = 1/(3*kB*T*VAVG)*trapz(JACF, dx=dt)/\
        units.length**3 * e**2 * units.velocity**2 * units.time
    return {jacf_out: {'t': TIME, 'JACF': JACF}, cond_out: COND}
