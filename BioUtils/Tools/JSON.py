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
Created on Jan 21, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import abc
import json
import gzip
try: import msgpack as mpk
except ImportError: mpk = None
from itertools import chain
from collections import Set
from functools import partial

class JSONable(object):
    __metaclass__ = abc.ABCMeta
    
    _separators = (',',':')
    _json_methods = [partial(json.dumps, separators=_separators), json.load] 
    methods = {'json' : _json_methods+[open], 
               'jgz' : _json_methods+[gzip.open],
               'json.gz' : _json_methods+[gzip.open],
               'mpk' : (mpk.dumps, mpk.load, open) if mpk is not None else (None, None) }
    default = 'jgz' if mpk is None else 'mpk'
    
    @abc.abstractmethod
    def jsonify(self): pass
    
    def _jsonify_attrs(self, attrs):
        return dict((attr, getattr(self, attr, None)) for attr in attrs)
    
    @classmethod
    def from_dict(cls, d): return cls(**d)
    
    def __str__(self): return str(self.jsonify())
    
    @staticmethod
    def _default(o):
        if isinstance(o, JSONable):
            return o.jsonify()
        if isinstance(o, Set):
            return tuple(o)
        raise TypeError(repr(o) + " is not JSON serializable")
    
    def save(self, filename, method=default):
        if method not in self.methods:
            raise ValueError('Wrong method name "%f". Available methods are: %s' 
                             % (method, self.methods.keys()))
        if mpk is None and method == 'mpk':
            print 'mpk is not available; using jgz instead. To use mpk install msgpack.'
            method = 'jgz'
        _dumps, _load, _open = self.methods[method]
        filename = filename+'.'+method
        out = _open(filename, 'wb')
        with out: out.write(_dumps(self, default=self._default))
        return filename
            
    @classmethod
    def load(cls, filename, method=None):
        if not method:
            method = filename[filename.index('.'):].strip('.')
        if method not in cls.methods:
            raise ValueError('Unsupported filetype "%s" (guessed by extension). Supported: %s' 
                             % (method, cls.methods.keys()))
        if mpk is None and method == 'mpk':
            print 'mpk is not available. To use mpk install msgpack.'
            return None
        _dumps, _load, _open = cls.methods[method]
        inp = _open(filename, 'rb')
        with inp: dct = _load(inp)
        assert isinstance(dct, dict), \
        '%s should contain a single top-level JSON object' % filename
        return cls.from_dict(dct)


class JSONslots(JSONable):
    def jsonify(self):
        return self._jsonify_attrs(chain.from_iterable(getattr(cls, '__slots__', []) 
                                                       for cls in type(self).__mro__))
            

#tests
if __name__ == '__main__':
    class Test(JSONslots):
        __slots__ = ['a', 'b']
        
        def __init__(self, a, b):
            self.a = a
            self.b = b
            
        @classmethod
        def from_dict(cls, d): return cls(**d)
        
    class Test1(Test):
        __slots__ = ['c', 'd']
        
        def __init__(self, a,b,c,d):
            super(Test1, self).__init__(a,b)
            self.c = c
            self.d = d
            
    t2 = Test1(1,2,3,4)
    print t2
            
    print
    t = Test((1,2,3,4,5), {'wert':1.35})
    fn = t.save('test')
    print fn
    with open(fn) as inp: print inp.read()
    print
    t1 = Test.load(fn)
    print t
    print t1
    
#    from BioUtils.Tools.tmpStorage import to_shelf, from_shelf
#    import cProfile as profile
#    import random, os, shutil
#    random.seed(42)
#    a = random.sample(xrange(1000000), 10000)
#    b = random.sample(xrange(1000000), 10000)
#    t = Test(a, b)
#    numtests = 100
#    profile.run('''for _i in xrange(numtests): 
#    t.save(fn)
#    Test.load(fn)
#    f = to_shelf(t)
#    from_shelf(f)''', 'json.profile')
#    f = to_shelf(t)
#    shutil.copy(f, os.path.basename(f))
    print 'Done'
    