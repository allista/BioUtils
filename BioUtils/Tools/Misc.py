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

def mktmp_name(suffix, register=True):
    fd, fn = tempfile.mkstemp(suffix, 'wb')
    if register: register_tmp_file(fn)
    os.close(fd)
    return fn
#end def

def safe_unlink(filename):
    try: os.unlink(filename)
    except: pass
