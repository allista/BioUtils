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
Created on Jan 20, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import sys, traceback
import psutil

from .Output import simple_timeit

class Pstats(object):
    def __init__(self, name=''):
        self.name = name
        self._cput = None
        self._p = None
        
    @staticmethod
    def _mgb(mem):
        dmem = vars(mem)
        return ' '.join('%s=%.3fMb' % (k, dmem[k]/1048576.) for k in dmem)

    def __enter__(self):
        self._p = psutil.Process()
        self._cput = self._p.cpu_times()
        print '[%d] S %s:\n%s\n' % (self._p.pid, self.name or self._p.name(), self._mgb(self._p.memory_info_ex()))
            
    def __exit__(self, *exc_info):
        cput1 = self._p.cpu_times()
        print '[%d] E %s:\n%s\n%s\n' % (self._p.pid, self.name or self._p.name(), self._mgb(self._p.memory_info_ex()),
                                        dict(user=cput1[0]-self._cput[0], sys=cput1[1]-self._cput[1]))
        
    def __call__(self, func):
        def ps_func(*args, **kwargs):
            old_name = self.name
            if not self.name: self.name = str(func)
            with self: ret = func(*args, **kwargs)
            self.name = old_name
            return ret
        copy_name(ps_func, func)
        return ps_func
    
class Timeit(object):
    def __init__(self, name=''):
        self.name = name
        
    def __call__(self, func):
        def t_func(*args, **kwargs):
            old_name = self.name
            if not self.name: self.name = str(func)
            with simple_timeit(self.name): 
                ret = func(*args, **kwargs)
            self.name = old_name
            return ret
        return t_func

def copy_name(func, obj):
    func.__name__ = getattr(obj, '__name__', obj.__class__.__name__)

def raise_tb(): raise RuntimeError(''.join(traceback.format_exception(*sys.exc_info())))

def raise_tb_on_error(func):
    def wrapper(*args, **kwargs):
        try: return func(*args, **kwargs)
        except: raise_tb()
    copy_name(wrapper, func)
    return wrapper