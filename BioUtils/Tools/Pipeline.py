'''
Created on Jan 30, 2016

@author: Allis Tauri <allista@gmail.com>
'''

from .Misc import safe_unlink
from .Output import isatty
import os, inspect

class PipelineNode(object):
    def __init__(self, func):
        self.outfile  = None
        self.inputs   = []
        self.required = []
        self.msg      = 'run %s' % func.__name__
        self.func     = func
        self.instance = None
        self.iclass   = None
        self.always_run = False
        self.argspec = inspect.getargspec(func)
        if self.argspec.keywords:
            print 'Warning: %s uses %s keywords argument. ' % func.__name__, self.argspec.keywords
        
    def nargs(self):
        nargs = len(self.argspec.args)
        if self.outfile: nargs -= 1
        if self.instance or self.iclass: nargs -= 1
        return nargs
        
    @classmethod
    def new(cls, required=[], inputs=[], msg=None, always_run=False):
        def factory(func):
            node = cls(func)
            node.inputs.extend(inputs)
            node.required.extend(required)
            node.always_run = always_run
            if msg: node.msg = msg
            return node
        return factory
    
    def copy(self):
        node = PipelineNode(self.func)
        node.msg = self.msg
        node.instance = self.instance
        node.iclass = self.iclass
        node.argspec = self.argspec
        node.always_run = self.always_run
        return node
    
    def plug(self, other):
        self.inputs.append(other)
        
    def unplug(self, other):
        try: self.inputs.remove(other)
        except: pass
        
    def replug(self, old, new):
        self.unplug(old)
        self.plug(new)
    
    @staticmethod
    def _isfile(filename):
        return (filename and 
                os.path.isfile(filename) and
                os.stat(filename).st_size > 0)
    
    @classmethod
    def _file_exists(cls, filename, msg):
        if cls._isfile(filename):
            print '\n%s exists and is up to date.' % filename
            if isatty:
                inp = raw_input('Are you sure you want to %s again? [y/n]: ' % msg)
                if inp.lower() != 'y': return True
            else:
                print 'To %s again, remove/rename the "%s" file.' % (msg, filename)
                return True
        return False
        
    def __get__(self, instance, owner):
        self.instance = instance
        self.iclass = owner
        return self
    
    def __call__(self, *inputs):
        inputs = self.inputs+list(inputs)
        nargs = self.nargs()
        if self.argspec.varargs: 
            if len(inputs) <= nargs:
                raise RuntimeError('PipelineNode %s needs at least %d input nodes, has %d' %
                                   (self.func.__name__, nargs+1, len(inputs)))
        elif len(inputs) != nargs:
            raise RuntimeError('PipelineNode %s needs exactly %d input nodes, has %d' %
                               (self.func.__name__, nargs, len(inputs)))
        args = [p.outfile for p in inputs]
        if args and self.outfile:
            print '\n'+', '.join(args)+' -> '+self.outfile
        if self.outfile: args = [self.outfile] + args
        force_run = False
        for inp in inputs:
            if self._isfile(inp.outfile):
                if self._isfile(self.outfile) \
                and os.stat(inp.outfile).st_mtime > os.stat(self.outfile).st_mtime:
                    force_run = True
            else:
                print '\nNo "%s" file found.' % inp.outfile
                if not inp(): return False
                force_run = True
        if not self.always_run and not force_run \
        and self._file_exists(self.outfile, self.msg): return True
        for r in self.required: r()
        if not self.func.__get__(self.instance, self.iclass)(*args):
            safe_unlink(self.outfile)
            return False
        return True
#end class