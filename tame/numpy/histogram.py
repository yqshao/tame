#!/usr/bin/env python
import numpy as np
from tame import FrameArray

class histogram(FrameArray):
    """histogram"""
    def __init__(self, var, *args, **kwargs):
        assert isinstance(var, FrameArray)
        self.var = var
        var.parent.derived.append(self)
        self.args = args
        self.kwargs = kwargs
        self.parent = var.parent
        self.update()

    def update(self):
        val = np.histogram(self.var.eval(), *self.args, **self.kwargs)
        self.val = val

    def eval(self):
        return self.val
