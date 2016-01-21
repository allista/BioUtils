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
Created on Mar 4, 2014

@author: Allis Tauri <allista@gmail.com>
'''

import sys
import traceback
from time import time
from datetime import timedelta
from contextlib import contextmanager
from multiprocessing.queues import Queue

isatty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()

class OutIntercepter(object):
    '''A file-like object which intercepts std-out/err'''
    def __init__(self):
        self._oldout = None
        self._olderr = None
    #end def
    
    def write(self, text): pass
    
    def flush(self): pass
    
    def __enter__(self):
        self._oldout = sys.stdout
        self._olderr = sys.stderr
        sys.stdout = sys.stdout = self
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        if _type is not None:
            print _value
            traceback.print_exception(_type, _value, _traceback, file=self._olderr)
        sys.stdout = self._oldout
        sys.stderr = self._olderr
        return True
    #end def
#end class


class OutQueue(OutIntercepter, Queue):
    '''A file-like object which puts text written into it 
    in a cross-process queue'''
    def __init__(self, maxsize=0):
        OutIntercepter.__init__(self)
        Queue.__init__(self, maxsize)
        self._id = None
    #end def

    def set_id(self, _id):
        self._id = _id
        return self
    #end def

    def put(self, obj, block=True, timeout=None):
        Queue.put(self, (self._id, obj), block=block, timeout=timeout)
        
    def write(self, text): self.put(text)
    
    def __enter__(self):
        OutIntercepter.__enter__(self)
        if self._id is None: self._id = -1
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        OutIntercepter.__exit__(self, _type, _value, _traceback)
        self._id = None
        return True
    #end def
#end class


@contextmanager
def user_message(msg, delimiter=' '):
    print msg+delimiter,
    sys.stdout.flush()
    try: yield
    finally:
        sys.stdout.write('Done\n')
        sys.stdout.flush()
        
        
class Progress(object):
    '''
    Simple tty progress indicator of the form
    > Message [ 56%]
    '''
    progress_len  = '      '
    progress_back = '\b'*len(progress_len)
    
    def __init__(self, msg, total, replace=True, is_index=True):
        self.start_msg = msg
        self.total     = total
        self.replace   = replace
        self.is_index  = is_index
        self._last_msg = ''
    
    def __enter__(self):
        if self.start_msg:
            if self.replace:
                print self.start_msg, self.progress_len if isatty else '',
            else: print self.start_msg
            sys.stdout.flush()
        return self
    
    def __exit__(self, *exc_info): 
        print ''
        sys.stdout.flush()
    
    def step(self, progress):
        if self.is_index: progress += 1
        if self.replace:
            msg = ('%s[%3.0f%%]' % 
                   (self.progress_back if isatty else ' ',
                    float(progress)/self.total*100.0))
            if msg != self._last_msg: sys.stdout.write(msg)
        else:
            msg = ('[%3.0f%%] %s of %s\n' % 
                   (float(progress)/self.total*100.0,
                    progress, self.total))
            sys.stdout.write(msg)
        self._last_msg = msg
        sys.stdout.flush()
#end class


class ProgressCounter(Progress):
    '''
    Simple tty progress indicator of the form
    > Message [ 56%]
    '''
    progress_len  = '      '
    progress_back = '\b'*len(progress_len)
    
    def __init__(self, msg, total, replace=True):
        assert(isinstance(total, int)), 'total should be integer'
        super(ProgressCounter, self).__init__(msg, total, replace, True)
        self._count = 0
    
    @property
    def percent(self):
        return float(self._count+1)/self.total*100.0
    
    def count(self):
        if self._count < self.total:
            self.step(self._count)
            self._count += 1
#end class

@contextmanager
def simple_timeit(name=''):
    try:
        t0 = time()
        print name+'...\n'
        yield 
    finally:
        print '%sElapsed time: %s\n' % (name+': ' if name else '', timedelta(seconds=time()-t0))
