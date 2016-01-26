#!/usr/bin/python
# coding=utf-8
#
# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# degen_primer is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
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