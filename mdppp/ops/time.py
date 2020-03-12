""" Time related functions"""
import numpy as np
from mdppp import FrameArray

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

def tcache(var, cache_size, **kwargs):
    """Cache of previous cache_size frames of data"""    
    return TCACHE(var, cache_size, **kwargs)

def tavg(var, drop_nan='all'):
    """Time Average Value"""
    return TAVG(var, drop_nan)

def tacf(var, cache_size, dropnan='partial'):
    var_cache = tcache(var, cache_size)
    acf = var_cache * var[np.newaxis, :]
    acf = np.mean(np.sum(acf, axis=2), axis=1)
    tacf = tavg(acf, dropnan)
    return tacf
    
def tmsd(var, cache_size, dropnan='partial'):
    var_cache = tcache(var, cache_size)
    dis = var_cache - var[np.newaxis, :]
    msd = np.mean(np.sum(dis**2, axis=2), axis=1)
    tmsd = tavg(msd, dropnan)
    return tmsd
