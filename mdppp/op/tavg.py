from mdppp import FrameArray

class TAVG(FrameArray):
    """
    """
    
    def __init__(self, var):
        assert type(var) == FrameArray
        FrameArray.__init__(self, var.parent)
        self.var = var
        self.cumsum = 0
        self.count = 0
        self.update()

    def update(self):
        self.cumsum += self.var.eval()
        self.count += 1

    def eval(self):
        return self.cumsum/self.count
    
def tavg(var):
    """Time Average Value

    Args:
        var: A FrameArray

    Returns:
        A TAVG Operation
    """
    return TAVG(var)
    
