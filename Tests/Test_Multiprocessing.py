#!/usr/bin/python
# coding=utf-8

'''
Created on Jan 25, 2016

@author: allis
'''

def setup(): 
    pass

def teardown(): 
    pass

def test():
    import multiprocessing as mp
    from BioUtils.Tools.Multiprocessing import MultiprocessingBase 
    abort_event = mp.Event()
    
    class Test(MultiprocessingBase):
        def __init__(self, abort_event):
            MultiprocessingBase.__init__(self, abort_event)
            self._a = 5
            
        def _method(self, i):
            return self._a+abs(i)
        
        def map_data(self, data):
            return self.parallelize_work(1, self._method, data)
        
        def map_functions(self, funcs, *args):
            return self.parallelize_functions(1, funcs, *args)
        
        def map_both(self, funcs, data):
            return self.parallelize_both(1, funcs, data)
    
    t = Test(abort_event)
    print t.map_data([-1, -2, -3, -4, -5, -6, -7, -8])
    print t.map_functions((abs, lambda x: x**2), -2)
    print t.map_both((abs, lambda x: x**2), (-2, 4))
    
if __name__ == '__main__':
    test()