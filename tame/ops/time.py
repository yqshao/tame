""" Time related functions"""
import numpy as np
from tame import FrameArray

class TCACHE(FrameArray):
    """Time     """
    def __init__(self, var, cache_size, default_val=np.nan):
        assert isinstance(var, FrameArray)
        self.var = var
        var.parent.derived.append(self)
        self.parent = var.parent       
        shape = var.eval().shape
        self.cache = np.full([cache_size, *shape], default_val)
        self.update()
    
    def update(self):
        new_val = np.expand_dims(self.var.eval(), 0)
        self.cache = np.concatenate(
            [new_val, self.cache], axis=0)[:-1]

    def eval(self):
        return self.cache

class TAVG(FrameArray):
    """Time avarage of FrameArray
    """
    
    def __init__(self, var, dropnan='all'):
        assert isinstance(var, FrameArray)
        var.parent.derived.append(self)
        self.var = var
        self.parent = var.parent
        self.dropnan = dropnan        
        self.cumsum = 0
        self.count = 0
        self.update()

    def update(self):
        new_val = self.var.eval()
        nan_val = np.isnan(new_val)
        if nan_val.any() and self.dropnan == 'all':
            return
        elif self.dropnan == 'partial':
            self.cumsum += np.nan_to_num(new_val)
            self.count += (1.0 - nan_val.astype(float))
            return
        self.cumsum += new_val
        self.count += 1

    def eval(self):
        return self.cumsum/self.count

class TSURVIVE(FrameArray):
    """Survival state FrameArray
    """

    def __init__(self, birth, death):
        assert isinstance(birth, FrameArray)
        assert isinstance(death, FrameArray)
        birth.parent.derived.append(self)
        self.birth = birth
        self.death = death
        self.parent = birth.parent
        self.val = birth.eval()

    def update(self):
        self.val = ((self.val | self.birth.eval()) & ~self.death.eval())

    def eval(self):
        return self.val

def tsurvive(birth, death):
    """Cache of previous cache_size frames of data

    Args:
        birth (FrameArray): boolean for birth
        death (FrameArray): boolean for death

    Returns:
        FrameArray: time-dependent survival state
    """
    return TSURVIVE(birth, death)

def tcache(var, cache_size, **kwargs):
    """Cache of previous cache_size frames of data

    Args:
        var (FrameArray): variable to be averaged
        cache_size (int): maximum $dt$ to calculate, in frames

    Returns:
        FrameArray: cached trajectory with shape [cache_size, input_shape]
    """
    return TCACHE(var, cache_size, **kwargs)

def tavg(var, dropnan='partial'):
    """Time average value:

    $$\langle \mathbf{var} \\rangle$$

    `tavg` ignores `nan` values when computing the average. When dropnan is
    'all', the entire `var` is discards if there is a nan inside. When dropnan
    is 'partial', the only `nan` entries are ignored.


    Args:
        var (FrameArray): variable to be averaged
        dropnan (str): 'all' or 'partial'

    Returns:
        FrameArray: time average of `var`
    """
    return TAVG(var, dropnan)

def tacf(var, cache_size, dropnan='partial'):
    """Time autocorrelation function:

    $$\langle \mathbf{var}_i(0) \cdot \mathbf{var}_i(dt) \\rangle_{i,t}$$

    `tacf` assumes that var has dimension `[n_particle, n_dim]`. The particle
    dimension will be averaged.

    Args:
        var (FrameArray): variable with shape `[n_particle, n_dim]`
        cache_size (int): maximum $dt$ to calculate, in frames
        dropnan (str): 'all' or 'partial'

    Returns:
        FrameArray: Autocorrelation function with shape `[cache_size]`
    """
    var_cache = tcache(var, cache_size)
    acf = var_cache * var[np.newaxis, :]
    acf = np.mean(np.sum(acf, axis=2), axis=1)
    tacf = tavg(acf, dropnan)
    return tacf
    
def tmsd(var, cache_size, dropnan='partial'):
    """Mean squred displacement function:

    $$\langle (\mathbf{var}_i(dt) - \mathbf{var}_i(0))^2 \\rangle_{i,t}$$

    `tmsd` assumes that var has dimension `[n_particle, n_dim]`. The particle
    dimension will be reduced by averaging.

    Args:
        var (FrameArray): variable with shape `[n_particle, n_dim]`
        cache_size (int): maximum $dt$ to calculate, in frames
        dropnan (str): 'all' or 'partial'

    Returns:
        FrameArray: Autocorrelation function with shape `[cache_size]`
    """
    var_cache = tcache(var, cache_size)
    dis = var_cache - var[np.newaxis, :]
    msd = np.mean(np.sum(dis**2, axis=2), axis=1)
    tmsd = tavg(msd, dropnan)
    return tmsd

def tpairsurvive(data, tag, n_cache, d_cache={}):
    """shorthand for the survival probability of pairs

    A pair is defined as two atoms of specific types (t1, t2) within certain
    distince criterion. Possible definitions of pair survival probabilities are:

    IMM (Impey, Madden, McDonald) is defined with one distance
    ---------------------------------------------------
    The pairs are defined as two ions that doees not depart for more than
    certain time takes a cutoff distance and a time (in ps) as argument.

    SSP (stable state picture) is defined with two distances
    --------------------------------------------------
    The correlation function is given by the probability of originating from a
    stable reactant and not reaching a stable product, takes two cutoff
    distances as arguments.

    The function returns a correlation function (frameArray with shape [dt, n1,
    n2]). at each time step, it gives a 0/1 state of whether a pair have been
    alive from dt ago.

    Args:
        data: FrameData
        tag: pair survival definition
        n_cache: number of frames to compute the correlation function
        d_cache(dict, optional): cached distances

    Returns:
        FrameArray: the survival correlation function
    """
    method = tag.split(':')[1]
    t1, t2 = map(int, tag.split(':')[0].split(','))
    rc = [*map(float, tag.split(':')[2].split(','))]
    coord = data['coord']
    cell = data['cell'][None, None,:]
    idx_i = np.where((data['elems']==t1).eval())[0]
    idx_j = np.where((data['elems']==t2).eval())[0]
    r_i = coord[idx_i]; r_j = coord[idx_j]
    if (t1, t2) not in d_cache:
        r_ij = r_i[:,None,:] - r_j[None, :, :]
        r_ij = r_ij - np.rint(r_ij/cell)*cell
        dist = np.sqrt(np.sum(r_ij**2, axis=2))
        d_cache[(t1, t2)] = dist
    birth = d_cache[(t1, t2)]<=rc[0]
    death = d_cache[(t1, t2)]>rc[int(method=='SSP')]
    if t1==t2:
        birth = birth * (1-np.eye(len(idx_i))[None, :, :])
    if method == 'IMM':
        n_death = int(rc[1]/dt)
        death = np.multiply.accumulate(
            tcache(death, n_death), axis=0)
        corrs = np.multiply.accumulate(
            tcache(tsurvive(birth, death), n_cache), axis=0)
    elif method == 'SSP':
        corr = np.multiply(
                tcache(birth, n_cache), np.multiply.accumulate(
                tcache(~death, n_cache), axis=0))#
    return corr
