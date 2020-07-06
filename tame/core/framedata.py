from .framearray import FrameArray

class FrameData():
    def __init__(self, generator):
        self.count = 0
        self.generator = generator                        
        data = next(self.generator)
        assert type(data) is dict, "Generator must return a dictionary"
        self.data = {k: FrameArray(self, val=v) for k, v in data.items()}
        self.derived = []

    def __repr__(self):
        repr = (f"<frameData at {self.count}: "+
                ", ".join(self.data.keys())+f"; {len(self.derived)} derived>")
        return repr
        
    def __getitem__(self, key):
        return self.data[key]

    def next_frame(self):
        data = next(self.generator)
        self.count += 1
        for k, v in data.items():
            self.data[k].val = v
        for child in self.derived:
            child.update()

    def run(self, steps=None):
        for i in range(steps):
            self.next_frame()
