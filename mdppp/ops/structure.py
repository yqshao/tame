"""Operations for structure manipulations"""

class UnwrappedCoord():
    def __init__(self, coord, cell):
        self.coord = coord
        self.cell = cell
        self.old_coord = coord.eval()
        self.val = self.old_coord
        
    def update():
        self.old_coord = self.val
        new_coord = self.coord.eval()
        self.val = old_coord + np.rint((new_coord-self.old_coord)/cell)*cell
    
    def eval():
        return self.val
        
def unwrap(coord, cell):
    return UnwrappedCoord(coord, cell)
