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
import tempfile
import signal
import traceback
import logging
from time import time, sleep
from datetime import timedelta
from contextlib import contextmanager
from datetime import datetime
from multiprocessing import Event

from BioUtils.Tools.tmpStorage import clean_tmp_files

re_type = type(re.compile(''))
isatty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
    
def retry(func, error_msg, num_retries):
    for i in xrange(num_retries):
        try:
            result = func()
            break
        except Exception as e:
            print(e)
            if i == num_retries-1:
                raise RuntimeError(error_msg)
    return result
#end def

@contextmanager
def user_message(msg, delimiter=' '):
    print msg+delimiter,
    sys.stdout.flush()
    try: yield
    finally:
        sys.stdout.write('Done\n')
        sys.stdout.flush()

def random_text(length):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _unused in xrange(length))

def mktmp_name(suffix):
    fd, fn = tempfile.mkstemp(suffix, 'wb')
    os.close(fd)
    return fn
#end def

def safe_unlink(filename):
    try: os.unlink(filename)
    except: pass


class MPMain(object):
    def __init__(self, pid):
        self.pid = pid
        self.abort_event = Event()
        signal.signal(signal.SIGINT,  self.sig_handler)
        signal.signal(signal.SIGTERM, self.sig_handler)
        signal.signal(signal.SIGQUIT, self.sig_handler)
        
    def sig_handler(self, signal, frame):
        if self.pid != os.getpid(): return
        print('\nAborting. This may take some time '
              'as not all operations could be stopped immediately.\n')
        self.abort_event.set(); sleep(0.1)
        clean_tmp_files()
    
    def _main(self): pass
    
    def __call__(self, *args, **kwargs):
        try: self._main()
        except:
            self.abort_event.set()
            traceback.print_exc()
            return 1
        return 0
#end class


class OutIntercepter(object):
    '''A file-like object which intercepts std-out/err'''
    def __init__(self):
        self._oldout = None
        self._olderr = None
    #end def
    
    def write(self, text): pass
    
    def flush(self): pass
    
    def __enter__(self):
        self._oldout = sys.stdout
        self._olderr = sys.stderr
        sys.stdout = sys.stdout = self
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        if _type is not None and not SystemExit:
            print(_value)
            traceback.print_exception(_type, _value, _traceback, file=self._olderr)
        sys.stdout = self._oldout
        sys.stderr = self._olderr
        return True
    #end def
#end class


class EchoLogger(OutIntercepter):
    '''Wrapper around logging module to capture stdout-err into a log file
    while still print it to std'''

    def __init__(self, name, level=logging.INFO):
        OutIntercepter.__init__(self)
        self._name      = name
        self._log       = name+'.log'
        self._level     = level
        self._logger    = logging.getLogger(name)
        self._handler   = logging.FileHandler(self._log, encoding='UTF-8')
        self._formatter = logging.Formatter('[%(asctime)s] %(message)s')
        #set log level
        self._handler.setLevel(self._level)
        self._logger.setLevel(self._level)
        #assemble pipeline
        self._handler.setFormatter(self._formatter)
        self._logger.addHandler(self._handler)
    #end def
        
    def __del__(self):
        self._handler.close()
        
        
    def __enter__(self):
        OutIntercepter.__enter__(self)
        self._logger.log(self._level, '=== START LOGGING ===')
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        if _type is not None and not SystemExit:
            print(_value)
            self._logger.error('Exception occured:', 
                               exc_info=(_type, _value, _traceback))
            traceback.print_exception(_type, _value, _traceback, file=self._olderr)
        sys.stdout = self._oldout
        sys.stderr = self._olderr
        self._logger.log(self._level, '=== END LOGGING ===')
        return True
    #end def
        
    def write(self, text):
        log_text = ' '.join(text.split())
        if log_text: self._logger.log(self._level, log_text, exc_info=False)
        self._oldout.write(text)
    #end def
    
    def flush(self): self._oldout.flush()
#end class


class Progress(object):
    '''
    Simple tty progress indicator of the form
    > Message [ 56%]
    '''
    progress_len  = '      '
    progress_back = '\b'*len(progress_len)
    
    def __init__(self, msg, total, replace=True, is_index=True):
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
    
    def step(self, progress):
        if self.is_index: progress += 1
        if self.replace:
            msg = ('%s[%3.0f%%]' % 
                   (self.progress_back if isatty else ' ',
                    float(progress)/self.total*100.0))
            if msg != self._last_msg: sys.stdout.write(msg)
        else:
            msg = ('[%3.0f%%] %s of %s\n' % 
                   (float(progress)/self.total*100.0,
                    progress, self.total))
            sys.stdout.write(msg)
        self._last_msg = msg
        sys.stdout.flush()
#end class


class ProgressCounter(Progress):
    '''
    Simple tty progress indicator of the form
    > Message [ 56%]
    '''
    progress_len  = '      '
    progress_back = '\b'*len(progress_len)
    
    def __init__(self, msg, total, replace=True):
        assert(isinstance(total, int)), 'total should be integer'
        super(ProgressCounter, self).__init__(msg, total, replace, True)
        self._count = 0
    
    @property
    def percent(self):
        return float(self._count+1)/self.total*100.0
    
    def count(self):
        if self._count < self.total:
            self.step(self._count)
            self._count += 1
#end class