
'''
Created on Feb 17, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from abc import ABCMeta

class DummyEvent(object):
    def is_set(self): return False
    def set(self): pass
#end class

def aborted(abort_event):
    try: return abort_event.is_set()
    except IOError: return True

class AbortableBase(object):
    '''Base class for all that need to be cleanly aborted through an event'''
    __metaclass__ = ABCMeta

    def __init__(self, abort_event):
        self._abort_event = abort_event
        self._aborted     = False
    #end def
    
    def aborted(self):
        if self._aborted: return True
        self._aborted = aborted(self._abort_event)
        return self._aborted
    #end def
#end class
        