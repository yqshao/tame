#!/usr/bin/env python3
import click
import tame
from tame.recipes.utils import load_traj_seg

@click.command(name='rdf',
               options_metavar='[options]',
               short_help='radial distribution function')
@load_traj_seg # general input handler
@click.option('--rmax', metavar='', default=30.0, show_default=True)
@click.option('--rbin', metavar='', default=0.1,  show_default=True)
@click.option('-t', '--tags', metavar='', default='1,1', show_default=True)
@click.option('--rdf-out',   metavar='', default='rdf', show_default=True)
def rdf_cmd(seg, dt, rmax, rbin, tags, rdf_out):
    """Computes the radial distribution function (RDF).

    See the documentation for detailed descriptions of the command:
    https://Teoroo-CMC.github.io/tame/latest/recipe/rdf
    """
    import numpy as np
    from tame import math as tm
    from tame.ops import tavg, build_nl
    from tame.io import load_traj
    # preparation
    tags = tags.split(' ')
    rgrid = np.arange(0, rmax, rbin)
    i_types = [int(t.split(',')[0]) for t in tags]
    j_types = [int(t.split(',')[1]) for t in tags]
    nl = build_nl(seg, i_types, j_types, rc=rmax)
    count, factor = {}, {}
    # setup counting FrameArrays
    for tag, i, j in zip(tags, i_types, j_types):
        _, _, dist = nl[i,j]
        icnt = np.count_nonzero(seg['elems'].val==i)
        jcnt = np.count_nonzero(seg['elems'].val==j)
        factor[tag] = icnt*jcnt
        count[tag] = tavg(tm.histogram(dist, bins=rgrid)[0])
    V = seg['cell'][0]*seg['cell'][1]*seg['cell'][2]
    Vavg = tavg(V)
    # run
    seg.run()
    # constants
    RDF = {'r': (rgrid[:-1]+rgrid[1:])/2}
    VAVG = Vavg.eval()
    vshell = 4*np.pi*RDF['r']**2*rbin
    # generating output
    for k, v in count.items():
        rho = factor[tag]/VAVG
        RDF[f'g_{{{k}}}(r)'] = v.eval()/rho/vshell
    return {rdf_out: RDF}
