# coding=utf-8

'''
Created on Jan 20, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os
import tempfile
from datetime import datetime
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

def run_cline(cline, *args, **kwargs):
    name = kwargs.pop('name', '')
    base = '%s%s-%s' % (name+'-' if name else '', 
                        os.path.basename(cline.program_name), 
                        datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    _raise = kwargs.pop('_raise', False)
    msg = kwargs.pop('_msg', None)
    err = kwargs.pop('stderr', base+'.err')
    out = kwargs.pop('stdout', base+'.out')
    print cline
    if isinstance(err, basestring): print 'stderr: %s' % err
    if isinstance(out, basestring): print 'stdout: %s' % out
    try: cline(stderr=err, stdout=out, *args, **kwargs)
    except Exception, e:
        if _raise: raise
        if msg: print msg
        print '%s\n' % str(e)
        return False
    if err and os.stat(err).st_size == 0:
        print 'Removing empty %s' % err 
        os.unlink(err)
    if out and os.stat(out).st_size == 0:
        print 'Removing empty %s' % err 
        os.unlink(out)
    return True

def mktmp_name(suffix='', register=True):
    fd, fn = tempfile.mkstemp(suffix, 'wb')
    if register: register_tmp_file(fn)
    os.close(fd)
    return fn

def safe_unlink(*filenames):
    for filename in filenames:
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
    
class IdNameMixin(object):
    @property
    def fullname(self):
        return ('%s (%s)' % (self.name, self.id) 
                if self.name != self.id else self.name) 
