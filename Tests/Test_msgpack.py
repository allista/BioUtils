'''
Created on Jan 19, 2016

@author: allis
'''

import msgpack
import pickle
from pydoc import locate
from DegenPrimer.SecStructures import Duplex

if __name__ == '__main__':
    d1 = Duplex('ATGCGCTA', 'ATGCGCTA', '1')
    d2 = Duplex('ATGCGCTA', 'ATGCGCTA', '2', revcomp=True)
    anns = (123, (d1, d2), {4:'assdf'})

    def packb_obj(obj):
        if hasattr(obj, '__dict__'):
            return dict(type=obj.__module__+'.'+type(obj).__name__, data=obj.__dict__)
        elif hasattr(obj, '__slots__'):
            return dict(type=obj.__module__+'.'+type(obj).__name__, data=dict((slot, getattr(obj, slot)) for slot in obj.__slots__))
        else: return str(obj)
            
    def unpackb_obj(obj):
        if 'type' in obj:
            t = locate(obj['type'])
            return t(**obj['data'])
        return obj
    
    b = msgpack.packb(anns, default=packb_obj)
    anns1 = msgpack.unpackb(b, use_list=False, object_hook=unpackb_obj)
    print anns1