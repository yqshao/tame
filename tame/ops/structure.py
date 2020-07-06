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
        
def unwrap(coord, cell):
    return UnwrappedCoord(coord, cell)
