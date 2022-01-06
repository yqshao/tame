#!/usr/bin/env python3
import click
import tame
from tame.recipes.utils import load_traj_seg

@click.command(name='persist',
               options_metavar='[options]',
               short_help='persistent time')
@load_traj_seg # general input handler
@click.option('--max-dt',     metavar='', default=30.0,          show_default=True)
@click.option('-t', '--tags', metavar='', default='1,1:SSP:3,4', show_default=True)
@click.option('--tcf-out',    metavar='', default='persist',     show_default=True)
def persist_cmd(seg, dt,
                tag, tcf_out, max_dt):
    """Computes the time correlation function

    See the documentation for detailed descriptions of the command:
    https://Teoroo-CMC.github.io/tame/latest/recipe/persist
    """
    import numpy as np
    import tame.math as tm
    from tame.ops import tavg, tpairsurvive
    from tame.io import load_traj

    n_cache = int(max_dt/dt)
    tags = tags.split(' ')

    corrs = {tag: # time correlation functions, one for each tag
             tavg(tm.mean(tm.sum(tpairsurvive(
                 data, tag, n_cache, d_cache), axis=2), axis=1))
             for tag in tags}
    seg.run()
    TIME = np.arange(0, n_cache)*dt
    CORRS = {'t': TIME}
    for k in keys:
        CORRS[k] = corrs[k].eval()
    return {tcf_out: CORRS}
