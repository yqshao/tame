"""Operations for structure manipulations"""
import numpy as np
from tame import FrameArray

class UnwrappedCoord(FrameArray):
    def __init__(self, coord, cell):
        assert isinstance(coord, FrameArray)
        coord.parent.derived.append(self)
        self.parent = coord.parent        
        self.coord = coord
        self.cell = cell
        self.val = coord.eval()
        
    def update(self):
        old_coord = self.val
        cell = self.cell.eval()
        new_coord = self.coord.eval()
        self.val = new_coord + np.rint((old_coord-new_coord)/cell)*cell
    
    def eval(self):
        return self.val

def _build_clist(ccount, cidx):
    natoms = cidx.shape[0]
    args  = np.lexsort([cidx[:,i] for i in range(3)])
    uc, uc0, ucinv, uccnt = np.unique(
        cidx[args], return_index=True, return_inverse=True, return_counts=True, axis=0)
    alist = np.full([ccount[0], ccount[1], ccount[2], uccnt.max()], -1)
    acid = np.zeros(natoms, int); acid[args] = np.arange(natoms)-uc0[ucinv]
    alist[cidx[:,0], cidx[:,1], cidx[:,2], acid] = np.arange(natoms)
    return uc, alist

class NeighbourList():
    def __init__(self, traj, types_i, types_j, rc):
        self.all_i = np.array([], np.int)
        self.all_j = np.array([], np.int)
        self.rc   = rc
        self.nl_i = FrameArray(traj, val=[])
        self.nl_j = FrameArray(traj, val=[])
        self.nl_r = FrameArray(traj, val=[])
        self.traj = traj
        for ti, tj in zip(types_i, types_j):
            if ti<=tj:
                ti, tj = tj, ti
            idx_i = np.where((traj['elems']==ti).eval())
            idx_j = np.where((traj['elems']==tj).eval())
            self.all_i = np.union1d(idx_i, self.all_i)
            self.all_j = np.union1d(idx_j, self.all_j)
        self.cell = traj['cell']
        self.coord_i = traj['coord'][self.all_i]
        self.coord_j = traj['coord'][self.all_j]
        self.all_i =  FrameArray(self.traj, val=self.all_i)
        self.all_j =  FrameArray(self.traj, val=self.all_j)
        traj.derived.append(self)
        self.update()

    def __getitem__(self, types):
        type1, type2 = types
        reverse = type1 > type2
        if reverse: type2, type1 = type1, type2
        idx_i, = np.where((self.traj['elems']==type1).eval())
        idx_j, = np.where((self.traj['elems']==type2).eval())
        i, j, r = self._get(idx_i, idx_j)
        if reverse: i, j = j, i
        return i, j, r

    def update(self):
        import itertools
        # builds the nl here
        def wrap(coord, cell):
            frac_coord = np.linalg.inv(cell.T)@coord.T
            frac_coord %= 1
            return (cell.T@frac_coord).T

        rc = self.rc
        cell = self.cell.val
        if np.ndim(cell) == 1: cell=np.diag(cell)
        coord_i = wrap(self.coord_i.val, cell)
        coord_j = wrap(self.coord_j.val, cell)
        all_i = self.all_i.val
        all_j = self.all_j.val
        n_repeat = rc * np.linalg.norm(np.linalg.inv(cell), axis=0)
        n_repeat = np.ceil(n_repeat)
        repeats = np.array([*itertools.product(*[np.arange(-n, n+1) for n in [1,1,1]])])
        coord_j_rep = coord_j[:,None,:] + repeats[None,:,:]@cell.T
        coord_j_rep = coord_j_rep.reshape([-1,3])
        idx_j_rep   = np.stack([np.arange(len(all_j))]*repeats.shape[0],
                               axis=1).reshape([-1])

        # create an array of cells
        cidx_i = np.floor(coord_i / rc).astype(int)
        cidx_j = np.floor(coord_j_rep / rc).astype(int)
        cimin = np.min(cidx_i, axis=0); cimax = np.max(cidx_i, axis=0)
        mask_j = (cidx_j<=(cimax+1)).all(axis=1) & (cidx_j>=(cimin-1)).all(axis=1)
        cidx_j = cidx_j[mask_j]; idx_j_rep = idx_j_rep[mask_j]; coord_j_rep = coord_j_rep[mask_j]
        cmin = np.min(np.concatenate([cidx_i, cidx_j], axis=0), axis=0)
        cidx_i -= cmin; cidx_j -= cmin
        ccount = np.max(np.concatenate([cidx_i, cidx_j], axis=0), axis=0)+1

        uc_i, alist_i = _build_clist(ccount, cidx_i)
        uc_j, alist_j = _build_clist(ccount, cidx_j)

        idx_i, idx_j, r = [], [], []
        for cidx in uc_i:
            alist_ci = alist_i[tuple(cidx)]
            alist_cj = alist_j[tuple(slice(cidx[i]-1,cidx[i]+2) for i in range(3))]
            alist_ci = alist_ci[alist_ci>-1]
            alist_cj = alist_cj[alist_cj>-1]
            coord_ci = coord_i[alist_ci]
            coord_cj = coord_j_rep[alist_cj]
            r_ij = np.linalg.norm(coord_ci[:,None,:] - coord_cj[None,:,:], axis=2)
            idx_i_rc, idx_j_rc = np.where(r_ij<rc)
            r_ij = r_ij[idx_i_rc, idx_j_rc]
            idx_i_rc = alist_ci[idx_i_rc]
            idx_j_rc = idx_j_rep[alist_cj[idx_j_rc]]
            idx_distinct, = np.where(all_i[idx_i_rc] != all_j[idx_j_rc])
            idx_i.append(idx_i_rc[idx_distinct])
            idx_j.append(idx_j_rc[idx_distinct])
            r.append(r_ij[idx_distinct])
        # # # updated the of i, j and r
        self.nl_i.val = np.concatenate(idx_i, axis=0)
        self.nl_j.val = np.concatenate(idx_j, axis=0)
        self.nl_r.val = np.concatenate(r, axis=0)

    def _get(self, idx_i, idx_j):
        """YS: for now NL is only accessed with static types of atoms,
        it might be desirable to allow nl.get for idx as FrameArrays.
        """
        bool_i = FrameArray(self.traj, val=np.in1d(self.all_i.val, idx_i))
        bool_j = FrameArray(self.traj, val=np.in1d(self.all_j.val, idx_j))
        mask = bool_i[self.nl_i] & bool_j[self.nl_j]
        return self.all_i[self.nl_i[mask]], self.all_j[self.nl_j[mask]], self.nl_r[mask]


def wrap(coord, cell):
    frac_coord = np.linalg.solve(cell.T, coord.T)
    frac_coord %= 1
    return (cell.T@frac_coord).T


def unwrap(coord, cell):
    return UnwrappedCoord(coord, cell)


def build_nl(traj, i_types, j_types, rc):
    return NeighbourList(traj, i_types, j_types, rc)
