"""Misc utils for recipes"""
import click

class seg_w_bar():
    def __init__(self, traj, nseg):
        self.traj = traj
        self.nseg = nseg
        self.iseg = 1
        self.derived = self.traj.derived

    def __getitem__(self, key):
        return self.traj[key]

    def run(self):
        with click.progressbar(range(self.nseg), label=f'seg{self.iseg}') as bar:
            for i in bar:
                self.traj.next_frame()
        self.iseg += 1
        self.traj.derived.clear()

def result_writer(results):
    import numpy as np
    # results should be a list of dicts
    for k, v in results[-1].items():
        if type(v) is dict:
            header = 'fields: ' + ' '.join(v.keys())
            if all([np.ndim(val)==0 for val in v.values()]):
                data = np.stack([[*r[k].values()] for r in results], axis=0)
                np.savetxt(f'{k}.dat', data, header=header)
            if all([np.ndim(val)==1 for val in v.values()]):
                data = np.stack([*v.values()], axis=1)
                np.savetxt(f'{k}_seg{len(results)}.dat', data, header=header)
        if np.ndim(v==0):
            # write all results to a same file if it's a scalar
            np.savetxt(f'{k}.dat', [r[k] for r in results])

def load_traj_seg(func):
    from tame.io import load_traj
    from functools import wraps

    @click.argument('trajs', metavar='trajs', nargs=-1, required=True)
    @click.option('-top',    metavar='', default=None, help="[default: None]")
    @click.option('-dt',     metavar='', default=1,    show_default=True)
    @click.option('-s', '--seg',    metavar='', default=5000, show_default=True)
    @click.option('-f', '--format', metavar='', default='auto', show_default=True)
    @wraps(func)
    def func_seg(trajs, top, format, seg, dt, **kwargs):
        traj =  load_traj(trajs, top, format)
        nseg = int(seg/dt)
        seg  = seg_w_bar(traj, nseg)
        results = []
        while True:
            try:
                result = func(seg, dt, **kwargs)
                results.append(result)
            except StopIteration:
                print('Reached end of traj.')
                break
            result_writer(results)
    return func_seg
