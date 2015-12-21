# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# BioUtils is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# indicator_gddccontrol is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on Jul 20, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import os, sys, re
import string
import random
from contextlib import contextmanager

from Bio import SeqIO

from DegenPrimer.MultiprocessingBase import MultiprocessingBase

re_type = type(re.compile(''))
isatty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()

class SeqLoader(MultiprocessingBase):
    def __init__(self, abort_event):
        super(SeqLoader, self).__init__(abort_event)
        
    @staticmethod
    def _load_file(filename, schema):
        if not os.path.isfile(filename): return None
        try: return SeqIO.read(filename, schema)
        except Exception, e:
            print 'Unable to parse %s' % filename
            print e
            return None
        
    def load_dir(self, dirname, schema, namefilter=None):
        if isinstance(namefilter, str):
            namefilter = re.compile(namefilter)
        if isinstance(namefilter, re_type):
            flt = lambda n: namefilter.match(n)
        if hasattr(namefilter, '__call__'):
            flt = namefilter
        else: flt = lambda n: True
        files = [os.path.join(dirname, f) 
                 for f in os.listdir(dirname) if flt(f)]
        return self.load(files, schema)
    
    def load_files(self, files, schema):
        return [r for r in self.parallelize_work(1, self._load_file, files, schema) if r]
    
def load_dir(abort_event, dirname, schema, namefilter=None):
    loader = SeqLoader(abort_event)
    return loader.load_dir(dirname, schema, namefilter)

def load_files(abort_event, filenames, schema):
    loader = SeqLoader(abort_event)
    return loader.load_files(filenames, schema)
    
class Progress(object):
    '''
    Simple tty progress indicator of the form
    > Message [ 56%]
    '''
    progress_len  = '      '
    progress_back = '\b'*len(progress_len)
    
    def __init__(self, msg, total):
        self.start_msg = msg
        self.total = total
    
    def __enter__(self):
        print self.start_msg, self.progress_len if isatty else '',
        sys.stdout.flush()
        return self
    
    def __exit__(self, *exc_info): pass
    
    def step(self, progress):
        sys.stdout.write('%s[%3.0f%%]' % 
                         (self.progress_back if isatty else ' ',
                          (float(progress)+1)/self.total*100.0))
        sys.stdout.flush()
#end class
    
@contextmanager
def user_message(msg):
    print msg+' ',
    sys.stdout.flush()
    try: yield
    finally:
        sys.stdout.write('Done\n')
        sys.stdout.flush()

def random_text(length):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _unused in xrange(length))
