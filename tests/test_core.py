import numpy as np

def test_framedata():
    from tame import FrameData
    from tame.ops import tavg
    
    data = np.array([1,2,3,4,5])
    iterator = map(lambda x: {'a':x}, iter(data))
    fdata = FrameData(iterator)

    a = fdata['a'] # A FrameArray object
    b = a * 2      # it should be able to do arithmetics
    c = np.sin(a)  # it should accept numpy ufuncs
    d = tavg(a)    # A FrameOp object
    # The above variables should be updated automatically
    for i in range(4):
        assert a.eval()==data[i]
        assert b.eval()==data[i]*2
        assert c.eval()==np.sin(data[i])
        assert d.eval()==np.mean(data[:i+1])
        fdata.next_frame()



    
