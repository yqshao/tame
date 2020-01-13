"""Computes autocorrelation function of a FrameArray"""
import numpy as np
from mdppp import FrameArray

class TACF(FrameArray):
    """Time Auto-Correlation Function Operation"""
    def __init__(self, var, cache_size):
        # Check inputs
        shape = var.eval().shape
        assert type(var) == FrameArray
        assert len(shape) == 2
        FrameArray.__init__(self, var.parent)
        # Initialize variables
        self.var = var
        self.cumsum = np.zeros(cache_size)
        self.count = np.zeros(cache_size)
        self.cache = np.full([cache_size, *shape], np.nan)        
        self.update()

    def update(self):
        new_val = np.expand_dims(self.var.eval(), 0)
        self.cache = np.concatenate(
            [new_val, self.cache], axis=0)[:-1]
        new_acf = np.mean(np.sum(self.cache * new_val, axis=2), axis=1)
        if not np.isnan(new_acf).any():
            self.cumsum += new_acf
            self.count += 1

    def eval(self):
        return self.cumsum/self.count    


def tacf(var, cache_size):
    """Time Autocorrelation Function 

    It is assumed that the varible is given in a 2d array.  The first
    dimension will be treated as equivalen particles, i.e. results
    will be automatically averaged. The second will be treated as the
    dimension of dimensions, where the 2nd order normalization is
    performed.

    For instance, when calculating the partial pressure
    autocorrelation function for the shear viscosity, ``var`` should
    have the shape [3, 1] (for three independent partial pressures)

    The functions returns a TACF operation which will be every new
    frame. Evaluation of the TACF operation will return a (reduced and
    normalized) autocorraltion function which the shape [cache_size]

    Args:
       var: A 2d (particle x dim) FrameArray.
       cache_size: Time interval (in frames).

    Returns:
        A TACF Operation
    """
    return TACF(var, cache_size)
