import numpy as np

symbols = ['X','H','HE','LI','BE','B','C','N','O','F','NE','NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE','I','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','ES','FM','MD','NO','LR','RF','DB','SG','BH','HS','MT','DS','RG','CN','NH','FL','MC','LV','TS','OG']

def _gromacs_trr_generator(trj, top, resort=False):
    import MDAnalysis as mda
    u=mda.Universe(top,trj)
    elems=np.array([symbols.index(t) for t in u.atoms.types])
    for ts in u.trajectory:
        data = {
            'coord': np.array(ts._pos, np.float64),
            'cell': ts.dimensions[0:3],
            'elems': elems}
        if ts.has_velocities:
            data['speed'] = np.array(ts._velocities, np.float64)
        yield data

def load_gromacs_trr(trj, top):
    from tame import FrameData
    return FrameData(_gromacs_trr_generator(trj, top))
