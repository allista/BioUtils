# coding=utf-8

'''
Uninterruptable MultiProcessing

Created on Mar 2, 2014

@author: Allis Tauri <allista@gmail.com>
'''

import signal 
import multiprocessing as mp
from multiprocessing.managers import BaseManager

from .Debug import Pstats #test


def ignore_interrupt():
    signal.signal(signal.SIGINT,  signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGQUIT, signal.SIG_IGN)
#end def


class UProcess(mp.Process):
    '''Process that ignores interrupt, term and quit signals'''
    def run(self):
        ignore_interrupt()
        mp.Process.run(self)
#end class


class UManager(BaseManager):
    '''Multiprocessing Manager that ignores interrupt, term and quit signals'''
    def start(self, initializer=ignore_interrupt, initargs=()):
        Pstats(self.__class__.__name__)(BaseManager.start)(self, initializer=initializer, initargs=initargs)#test
#end class


class AutoProxyMeta(type):
    '''Metaclass that replicates multiprocessing.managers.MakeProxyType
    functionality, but allows proxy classes that use it to be pickable'''
    def __new__(cls, name, bases, attrs):
        dic = {}
        for meth in attrs.get('_exposed_', ()):
            exec '''def %s(self, *args, **kwds):
            return self._callmethod(%r, args, kwds)''' % (meth, meth) in dic
        dic.update(attrs)
        return super(AutoProxyMeta, cls).__new__(cls, name, bases, dic)
#end metaclass


class _FuncManager(UManager):
    _functions = set()
    
    @classmethod
    def functions(self): return self._functions
    
    @classmethod
    def _register_function(cls, name):
        #make a copy of parents' _functions
        if '_functions' not in cls.__dict__:
            cls._functions = cls._functions.copy()
        #create new method
        def temp(self, *args, **kwargs):
            func = getattr(cls, '_'+name, lambda x: None)
            result = func(self, *args, **kwargs)
            return result._getvalue()
        temp.__name__ = name
        setattr(cls, name, temp)
        cls._functions.add(name)
    #end def
    
    @classmethod
    def add_functions(cls, funcs, names=None):
        for i, func in enumerate(funcs):
            try: name = names[i]
            except: name = (func.__name__).lstrip('_')
            cls.register('_'+name, func)
            cls._register_function(name)
    #end def
#end class

def FuncManager(typeid, funcs, names=None):
    '''Return a subclass of UManager with funcs registered by names (if given)
    or by (func.__name__).lstrip('_') for func in funcs'''
    manager = type(typeid, (_FuncManager,), {})
    manager.add_functions(funcs, names)
    return manager
#end def


def at_manager(Manager, func_name):
    def temp(*args, **kwargs):
        mgr = Manager(); mgr.start()
        func = getattr(mgr, func_name)
        result = func(*args, **kwargs)
        mgr.shutdown()
        return result
    temp.__name__ = func_name
    return temp
#end def
