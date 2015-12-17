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

import sys
import string
import random
from contextlib import contextmanager
import multiprocessing

isatty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
ncpu = multiprocessing.cpu_count()

class Progress(object):

    progress_len  = '      '
    progress_back = '\b'*len(progress_len)
    
    @classmethod
    def step(cls, progress, total):
        sys.stdout.write('%s[%3.0f%%]' % 
                         (cls.progress_back if isatty else ' ',
                          (float(progress)+1)/total*100.0))
        sys.stdout.flush()
    
    @classmethod
    def start(cls, msg):
        print msg, cls.progress_len if isatty else '',
        sys.stdout.flush()
    
@contextmanager
def user_message(msg):
    print msg,
    sys.stdout.flush()
    try: yield
    finally:
        sys.stdout.write(' Done\n')
        sys.stdout.flush()

def random_text(length):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _unused in xrange(length))
