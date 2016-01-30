'''
Created on Jan 30, 2016

@author: Allis Tauri <allista@gmail.com>
'''

class Lineage(tuple):
    '''Handles lineage information'''
    
    _prokaryotes = ('archaea', 'bacteria')
    known_taxons = set(_prokaryotes)
    
    @classmethod
    def _register(cls, lineage):
        for t in lineage:
            cls.known_taxons.add(t)
    
    def __new__(cls, line, delimiter=';'):
        line = line.strip(' \n\r')
        taxons = [t.strip().lower() for t in line.split(delimiter)]
        if taxons and taxons[0] in cls._prokaryotes:
            taxons.insert(0, 'prokaryotes')
        return tuple.__new__(cls, taxons)
    
    def __init__(self, line, delimiter=';'):
        self.str = ':'.join(self)
        self._register(self)
        
    def __str__(self): return self.str
    
    def __eq__(self, other):
        return self.str == other.str
    
    def __sub__(self, other):
        return Lineage(self.str.replace(other.str, '', 1).strip(':'), delimiter=':')
    
    @property
    def last(self):
        return self[-1] if self else None
    
    def includes(self, other):
        return other.str.startswith(self.str)
    
    @classmethod
    def from_iter(cls, taxons):
        if not taxons: return cls('')
        return cls(';'.join(taxons))
    
    @classmethod
    def common(cls, lineages, delimiter=':'):
        if not lineages: return None
        i = 0
        common = []
        lineages = list(lineages)
        while True:
            clades = set()
            for l in lineages:
                if isinstance(l, str):
                    l = Lineage(l, delimiter)
                if not isinstance(l, Lineage):
                    raise ValueError('Unsupported lineage type: %s' % type(l))
                if len(l) <= i:
                    clades.clear() 
                    break
                if len(clades) > 0:
                    if clades.add(l[i]): break
                else: clades.add(l[i])
            if len(clades) != 1: break
            common.append(clades.pop())
            i += 1
        return cls(';'.join(common))
    
    @classmethod
    def sameall(cls, lineages):
        if not lineages: return False
        l0 = lineages[0]
        for l in lineages[1:]:
            if l0 != l: return False
        return True
    
    @classmethod
    def samelast(cls, lineages):
        if not lineages: return False
        l0 = lineages[0].last
        for l in lineages[1:]:
            if l0 != l.last: return False
        return True
#end class