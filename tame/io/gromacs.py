import MDAnalysis as mda
import numpy as np

symbols = ['X','H','HE','LI','BE','B','C','N','O','F','NE','NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE','I','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','ES','FM','MD','NO','LR','RF','DB','SG','BH','HS','MT','DS','RG','CN','NH','FL','MC','LV','TS','OG']

def _gromacs_trr_generator(trj, top, resort=False):
    u=mda.Universe(top,trj)
    for ts in u.trajectory:
        coord_ts=ts._pos
        speed_ts=ts._velocities
        elem_ts=[symbols.index(t) for t in u.atoms.types]
        dim=ts.dimensions
        cell=dim[0:3]
        coord, speed, elems = [], [], []
        for n in range(len(elem_ts)):
            coord.append(coord_ts[n])
            speed.append(speed_ts[n])
            elems.append(elem_ts[n])
        data = {
            'coord': np.array(coord, np.float64),
            'speed': np.array(speed, np.float64),
            'elems': np.array(elems),
            'cell': cell}
        yield data

def load_gromacs_trr(trj, top):
    from tame import FrameData
    return FrameData(_gromacs_trr_generator(trj, top))
