import numpy as np

def test_framedata():
    from mdppp import FrameData
    from mdppp.ops import CumMean
    
    data = np.array([1,2,3,4,5])
    iterator = iter(data)
    fdata = FrameData({'a': iterator})

    a = fdata['a'] # A FrameArray object
    b = a * 2      # it should be able to do arithmetics
    c = np.sin(a)  # it should accept numpy ufuncs
    d = CumMean(a) # A FrameOp object
    # The above variables should be updated automatically
    for i in range(5):
        assert a.eval()==data[i]
        assert b.eval()==data[i]*2
        assert c.eval()==np.sin(data[i])
        assert d.eval()==np.mean(data[:i])
        fdata.next_frame()



    
