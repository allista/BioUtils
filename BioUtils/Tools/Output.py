# coding=utf-8

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
        self._out = None
        self._err = None
    #end def
    
    def write(self, text): pass
    
    def flush(self): pass
    
    def __enter__(self):
        self._out = sys.stdout
        self._err = sys.stderr
        sys.stdout = sys.stdout = self
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        if _type is SystemExit: pass
        elif _type is not None:
            print _value
            traceback.print_exception(_type, _value, _traceback, file=self._err)
        sys.stdout = self._out
        sys.stderr = self._err
        return True
    #end def
#end class

class OutQueue(OutIntercepter, Queue):
    '''A file-like object which puts text written into it 
    in a cross-process queue'''
    def __init__(self, maxsize=0):
        OutIntercepter.__init__(self)
        Queue.__init__(self, maxsize)

    def write(self, text): self.put(text)
#end class

class OutQueueWithID(OutIntercepter, Queue):
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
    
    def __init__(self, msg, total=0, replace=True, is_index=True):
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
    
    def _progress_string(self, progress):
        return ('[%3.0f%%]' % (float(progress)/self.total*100.0) 
                if self.total > 0 else '[%s]' % progress)
    
    def step(self, progress):
        if self.is_index: progress += 1
        if self.replace:
            msg = ('%s%s' %
                   (self.progress_back if isatty else ' ',
                    self._progress_string(progress)))
            if msg != self._last_msg: sys.stdout.write(msg)
        else:
            msg = '%s' % self._progress_string(progress)
            if self.total > 0: msg += ' %s of %s' % (progress, self.total)
            sys.stdout.write(msg+'\n')
        self._last_msg = msg
        sys.stdout.flush()
#end class


class ProgressCounter(Progress):
    '''
    Simple tty progress indicator of the form
    > Message [ 56%]
    '''
    def __init__(self, msg, total=0, replace=True):
        assert(isinstance(total, int)), 'total should be integer'
        super(ProgressCounter, self).__init__(msg, total, replace, True)
        self._count = 0
    
    @property
    def percent(self):
        return (float(self._count+1)/self.total*100.0 
                if self.total > 0 else self._count+1)
    
    def count(self):
        if self.total == 0 or self._count < self.total:
            self.step(self._count)
            self._count += 1
#end class

class simple_timeit(object):
    def __init__(self, name=''):
        self.name = name
        self.timedelta = 0
        self.seconds = 0
        self._t0 = 0
        
    def __enter__(self):
        self._t0 = time()
        print '%s...\n' % self.name
        sys.stdout.flush()

    def __exit__(self, *exc_info):
        self.timedelta = timedelta(seconds=time()-self._t0)
        self.seconds = self.timedelta.total_seconds()
        print '%sElapsed time: %s\n' % (self.name+': ' if self.name else '', self.timedelta)
        sys.stdout.flush()