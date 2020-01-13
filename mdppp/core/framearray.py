import numpy as np
import mdppp

class FrameArray():
    """
    """
    
    def __init__(self, parent, val=None,
                 fn=None, fn_args=None, fn_kwargs=None):
        self.parent = parent
        self.val = val
        self.fn = fn
        self.fn_args = fn_args
        self.fn_kwargs = fn_kwargs
        if val is None:
            self.parent.derived.append(self)

    def __repr__(self):
        return "<FrameArray>"

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Given a function accepting some arguments, returning a updating
        function that the parent could call to update the derived
        varible.
        """
        new_child = FrameArray(self.parent,
                               fn=ufunc, fn_args=inputs, fn_kwargs=kwargs)
        new_child.update()
        return new_child
    
    def __getitem__(self, key):
        return self.__array_ufunc__(
            (lambda x, key: x[key]), None, self, key)

    def __lt__(self, other):
        return np.less(self, other)

    def __le__(self, other):
        return np.less_equal(self, other)

    def __eq__(self, other):
        return np.equal(self, other)

    def __ne__(self, other):
        return np.not_equal(self, other)

    def __gt__(self, other):
        return np.greater(self, other)
    
    def __ge__(self, other):
        return np.greater_equal(self, other)

    def __add__(self, other):
        return np.add(self, other)

    def __sub__(self, other):
        return np.subtract(self, other)

    def __mul__(self, other):
        return np.multiply(self, other)

    def __truediv__(self, other):
        return np.true_divide(self, other)

    def __floordiv__(self, other):
        return np.floor_divide(self, other)

    def __mod__(self, other):
        return np.mod(self, other)

    def __pow__(self, other):
        return np.power(self, other)

    def sum(self, **kwargs):
        return self.__array_ufunc__(np.sum, None, (self), **kwargs)

    def mean(self, **kwargs):
        return self.__array_ufunc__(np.mean, None, (self), **kwargs)

    def update(self):
        args = [arg.eval() if isinstance(arg, FrameArray)
                else arg for arg in self.fn_args]
        kwargs = {k: arg.eval() if isinstance(arg, FrameArray)
                   else arg for k, arg in self.fn_kwargs.items()}
        self.val = self.fn(*args, **kwargs)
        
    def eval(self):
        """Return the value of array"""
        return self.val
