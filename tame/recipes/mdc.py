#!/usr/bin/env python3
import click
import tame
import numpy as np
from tame.fit import linear_fit
from tame.ops import tmsd, tpairsurvive, tavg, tcache, unwrap
from tame.units import units, kB, e
from tame.recipes.utils import load_traj_seg # general innput handler

def make_mdc(data, tag, n_cache):
    all_dr = {}
    coord = data['coord']
    def get_dr(t): # dr is reconstructed at each call to make_corrs
        if t not in all_dr.keys():
            r_t = coord[data['elems'].eval()==t]
            all_dr[t] = tcache(r_t, n_cache) - r_t[None, :, :]
        return all_dr[t]

    if ',' not in tag: # self_diffusion, return MSD
        corr = np.mean(np.sum(get_dr(int(tag))**2, axis=2), axis=1)
        return {tag:corr}
    elif ':' not in tag: # distinct diffusion, return MCD
        t1, t2 = map(int, tag.split(','))
        dr_i = get_dr(t1); dr_j = get_dr(t2)
        if t1==t2: # same species
            dr_j = np.sum(dr_j, axis=1, keepdims=True) - dr_i
            corr = np.mean(np.sum(dr_i*dr_j, axis=2), axis=1)
        else: # different species
            corr = np.sum(
                np.mean(dr_i,axis=1)*np.sum(dr_j,axis=1), axis=1)
        return {tag:corr}
    else:
        t1, t2 = map(int, tag.split(':')[0].split(','))
        dr_i = get_dr(t1); dr_j = get_dr(t2)
        alive = tpairsurvive(data, tag, n_cache, d_cache)
        dead = -(alive-1)-np.eye(alive.val.shape[1])[None,:,:]*int(t1==t2)
        dr_j_in  = np.sum(dr_j[:,None,:,:]*alive[:,:,:,None],axis=2)
        dr_j_out = np.sum(dr_j[:,None,:,:]*dead[:,:,:,None],axis=2)
        cnt = np.mean(np.sum(alive, axis=2), axis=1)
        corr_in = np.mean(np.sum(dr_i*dr_j_in, axis=2), axis=1)
        corr_out = np.mean(np.sum(dr_i*dr_j_out, axis=2), axis=1)
        return {f'{tag}:cnt': tavg(cnt),
                f'{tag}:in':  tavg(corr_in),
                f'{tag}:out': tavg(corr_out)}

@click.command(name='mdc',
               options_metavar='[options]',
               short_help='mean displacement correlation')
@load_traj_seg # general input handler
@click.option('--max-dt',    metavar='', default=30.0,   show_default=True)
@click.option('-t', '--tag', metavar='', default='3',    show_default=True)
@click.option('--mdc-out',   metavar='', default='mdc',  show_default=True)
@click.option('--fit-min',   metavar='', default=5.0,    show_default=True)
@click.option('--fit-max',   metavar='', default=20.0,   show_default=True)
@click.option('--diff-out',  metavar='', default='diff', show_default=True)
def mdc_cmd(seg, dt, #
            tag, max_dt, fit_min, fit_max, mdc_out, diff_out):
    """Commputing the self/distinct diffusion coefficients using mean
    displacement correlations (MDC).

    See the documentation for detailed descriptions of the command:
    https://Teoroo-CMC.github.io/tame/latest/recipe/mdc
    """
    n_cache = int(max_dt/dt)
    coord = unwrap(seg['coord'], seg['cell'])

    corrs = {}
    for t in tag.split(' '):
        corrs.update(make_mdc(seg, t, n_cache))
    keys = list(corrs.keys())
    labels = [f'D_{k}' for k in list(corrs.keys())]

    seg.run()
    TIME = np.arange(0, n_cache)*dt
    CORRS = {'t': TIME}
    diffs = {}
    for k in keys:
        CORRS[k] = corrs[k].eval()
        diffs[k] = 1/6*linear_fit(TIME, CORRS[k], fit_min, fit_max)[0]
    return {mdc_out: CORRS, diff_out: diffs}
