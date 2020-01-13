"""Computes autocorrelation function of a FrameArray"""

class TMSD(FrameArray):
    """Time Mean Squared Displacement Operation"""

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
        new_val = np.expand_dims(val.eval(), 0)
        self.cache = np.concatenate(
            [self.cache, new_val], axis=0)[1:]
        new_msd = self.cache * new_val
        # self.cumsum += 
        # self.count += 

    def eval(self):
        return self.cumsum/self.count    


def tmsd(var, cache_size):
    """Time Mean Squared Displacement 

    It is assumed that the varible is given in a 2d array.  The first
    dimension will be treated as equivalen particles, i.e. results
    will be automatically averaged. The second will be treated as the
    dimension of dimensions, where the 2nd order normalization is
    performed.

    The functions returns a TMSD operation which will be every new
    frame. Evaluation of the TMSD operation will return a (reduced and
    normalized) autocorraltion function which the shape [cache_size]

    Args:
       var: A 2d (particle x dim) FrameArray.
       cache_size: Time interval (in frames).

    Returns:
        A TMSD Operation
    """
    return TMSD(var, cache_size)
