
'''
Created on Jan 21, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import abc
import json
import gzip
try: import msgpack as mpk
except ImportError: mpk = None

from collections import Set
from functools import partial

from pydoc import locate

    
class JSONable(object):
    """Abstract base class for custom classes than need to be serializable to JSON"""
    __metaclass__ = abc.ABCMeta

    _separators = (',',':')
    _json_methods = [partial(json.dumps, separators=_separators), json.load] 
    methods = {'json' : _json_methods+[open], 
               'jgz' : _json_methods+[gzip.open],
               'json.gz' : _json_methods+[gzip.open],
               'mpk' : (mpk.dumps, mpk.load, open) if mpk is not None else (None, None) }
    default_method = 'jgz'# if mpk is None else 'mpk'
    
    @classmethod
    def qualname(cls): 
        """Return a fully qualified name of a class.
        In python >= 3.3 it's a true qualified name; 
        In python < 3.3 only the top-level classes are supported"""
        return getattr(cls, '__qualname__', 
                   '%s.%s' % (cls.__module__,
                              cls.__name__))
        
    _cls_key = '__cls'
    _data_key = '__data'
    
    @classmethod
    def _pack(cls, name, obj):
        return {cls._cls_key: name, cls._data_key: obj}
        
    @abc.abstractmethod
    def jsonify(self): pass 
    
    @classmethod
    @abc.abstractmethod
    def from_dict(cls, d): pass
    
    def __str__(self): 
        return json.dumps(self, indent=4, default=self._jsonify_obj)
    
    @classmethod
    def _jsonify_obj(cls, o):
        if isinstance(o, JSONable): return cls._pack(o.qualname(), o.jsonify())
        raise TypeError('%s is not JSON serializable' % type(o))
    
    def save(self, filename, method=None):
        '''
        @filename: may be a basename or a full name with extension; in the latter
        case the extension will be used to guess method, if it is None
        @method: should be one of the keys of JSONable.methods or None
        '''
        #try to guess method from filename; or use default method
        if not method:
            try: dot = filename.index('.')
            except ValueError: dot = len(filename)
            method = filename[dot:].strip('.')
            if method not in self.methods: method = self.default_method
        #check provided method
        if method not in self.methods:
            raise ValueError('Wrong method name "%s". Available methods are: %s' 
                             % (method, self.methods.keys()))
        #check for msgpack
        if mpk is None and method == 'mpk':
            print 'mpk is not available; using jgz instead. To use mpk install msgpack.'
            method = 'jgz'
        #save and return a possibly new filename
        _dumps, _load, _open = self.methods[method]
        ext = '.'+method
        if not filename.endswith(ext): filename += ext
        out = _open(filename, 'wb')
        with out: out.write(_dumps(self, default=self._jsonify_obj))
        return filename
            
    @classmethod
    def load(cls, filename, method=None):
        '''
        @filename: a path to a file with the saved object; its extension will be 
        used to guess method if is None 
        @method: should be one of the keys of JSONable.methods or None
        '''
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
        assert dct.get(cls._cls_key) == cls.qualname(), \
        'Type mismatch: trying to load %s, but %s was saved' % (cls.qualname(), dct.get(cls._cls_key))
        return cls.from_dict(dct.get(cls._data_key, {}))


class JSONattrs(JSONable):
    """JSONable class that saves/loads all its attributes 
    (both from __dict__ and __slots__). If an attribute is a JSONable, it is
    recreated from its class on load. 
    """

    @classmethod
    def _pack_obj(cls, obj):
        if isinstance(obj, Set):
            return cls._pack(type(obj).__name__, tuple(obj))
        if isinstance(obj, tuple):
            return cls._pack('tuple', obj)
        return obj
    
    @classmethod
    def _unpack_obj(cls, obj):
        if isinstance(obj, list):
            return [cls._unpack_obj(it) for it in obj]
        if isinstance(obj, dict):
            T = obj.get(cls._cls_key)
            if not T: return obj
            T = locate(T)
            if not T: return obj
            data = obj.get(cls._data_key)
            if data is None: return T()
            if issubclass(T, JSONable):
                return T.from_dict(data)
            return T(cls._unpack_obj(data))
        return obj
    
    def jsonify(self):
        d = dict()
        for cls in type(self).__mro__:
            for attr_name in getattr(cls, '__slots__', []): 
                d[attr_name] = self._pack_obj(getattr(self, attr_name, None))
        if hasattr(self, '__dict__'):
            for attr_name, attr in self.__dict__.iteritems():
                d[attr_name] = self._pack_obj(attr)
        return d
    
    @classmethod
    def _default(cls):
        raise NotImplementedError(('In %s either _default() classmethod should be implemented '
                                   'or __init__ should be callable without arguments.') % cls.__name__) 
    
    @classmethod
    def default(cls): 
        try: return cls()
        except TypeError:
            return cls._default()
    
    @classmethod
    def from_dict(cls, d):
        o = cls.default() 
        for attr_name in d:
            setattr(o, attr_name, cls._unpack_obj(d[attr_name]))
        return o
            

#tests
if __name__ == '__main__':
    import binascii
    class Test(JSONattrs):
        __slots__ = ['a', 'b']
        
        def __init__(self, a, b):
            self.a = a
            self.b = b
            
        @classmethod
        def _default(cls): return cls(0,0)
            
    class Test1(Test):
        __slots__ = ['c', 'd']
        
        def __init__(self, a,b,c,d):
            super(Test1, self).__init__(a,b)
            self.c = c
            self.d = d
            
        @classmethod
        def _default(cls): return cls(0,0,0,0)
            
    class Test2(Test1):
        def __init__(self):
            super(Test2, self).__init__(1,2,3,4)
            self._t1 = (1,2,3,4,5)
            self._t2 = {'wert':1.35}
            
    class Test3(JSONattrs):
        def __init__(self, **kwargs):
            self._A = Test(kwargs.get('a'),kwargs.get('b'))
            self._B = Test1(kwargs.get('b'),kwargs.get('a'), 5, 9)
            self._C = Test2()
            
    t = Test3(a=3,b=7)
    fn = t.save('test')
    print fn
    with open(fn) as inp: print binascii.b2a_hex(inp.read())
    print
    t1 = Test3.load(fn)
    print t
    print
    print t1
    assert str(t) == str(t1)
    
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
    