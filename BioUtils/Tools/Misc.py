# coding=utf-8

'''
Created on Jan 20, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os
import tempfile
from .tmpStorage import register_tmp_file
    
def retry(func, error_msg, num_retries):
    for i in xrange(num_retries):
        try:
            result = func()
            break
        except Exception as e:
            print e
            if i == num_retries-1:
                raise RuntimeError(error_msg)
    return result
#end def

def mktmp_name(suffix='', register=True):
    fd, fn = tempfile.mkstemp(suffix, 'wb')
    if register: register_tmp_file(fn)
    os.close(fd)
    return fn
#end def

def safe_unlink(filename):
    try: os.unlink(filename)
    except: pass

class ListDB(object):
    def __init__(self):
        self._db = {}
    
    def get(self, key, default=None):
        return self._db.get(key, default)
    
    def __getitem__(self, key): 
        return self._db[key]
    
    def __setitem__(self, key, value):
        l = self._db.get(key, None)
        if l is None: self._db[key] = [value]
        else: l.append(value)
        
    def __iter__(self): return self._db.__iter__()
        