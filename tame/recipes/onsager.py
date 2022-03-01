#!/usr/bin/env python
import click
from tame.recipes.utils import load_traj_seg

@click.command(name='onsager',
               options_metavar='[options]',
               short_help='onsager coefficients')
@load_traj_seg # general input handler
@click.option('--max-dt',     metavar='', default=20.0,   show_default=True)
@click.option('--rcom',       metavar='', default=None,   show_default=True)
@click.option('-t', '--tags', metavar='', default='3',    show_default=True)
@click.option('--corr-out',   metavar='', default='corr', show_default=True)
def onsager_cmd(seg, dt, rcom, tags, max_dt, corr_out):
    """Computing the onsager coefficients with the mean displacement correlation
    between species.

    See the documentation for detailed descriptions of the command:
    https://Teoroo-CMC.github.io/tame/latest/recipe/onsager

    """
    import numpy as np
    from tame.fit import linear_fit
    from tame.ops import tavg, tcache, unwrap
    from tame.io import load_traj
    # Processing input
    n_cache = int(max_dt/dt)
    tags = tags.split(' ')
    i_types = [int(tag.split(',')[0]) for tag in tags]
    j_types = [int(tag.split(',')[1]) for tag in tags]
    dP = {} # cache the displacement for each type of interest

    # Build correlations
    coord = unwrap(seg['coord'], seg['cell'])
    if rcom is not None:
        am = {int(k):float(v) for k,v in map(lambda x:x.split(':'), rcom.split(','))}
        masses = np.array([am[e] for e in seg['elems'].eval()])
        com = np.sum(coord*masses[:,None], axis=0, keepdims=True)/masses.sum()
        coord = coord-com

    for t in np.union1d(i_types, j_types):
        t_idx = (seg['elems'] == t)
        P_t = np.sum(coord[t_idx], axis=0)
        dP[t] = tcache(P_t, n_cache) - P_t[None, :]

    corrs = {} # build the correlation functions here
    for tag, i, j in zip(tags, i_types, j_types):
        corrs[tag] = tavg(np.sum(dP[i]*dP[j], axis=1), dropnan='partial')

    # Run
    seg.run()

    # Generating output
    TIME = np.arange(0, n_cache)*dt
    CORRS = {'t': TIME}
    for tag in tags: # correlation functions are normalized wrt the volume!
        CORRS[tag] = corrs[tag].eval()

    return {corr_out: CORRS}
