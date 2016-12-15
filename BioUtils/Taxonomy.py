'''
Created on Jan 30, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import re

from SeqUtils import SeqView, pretty_rec_name

class Lineage(tuple):
    '''Handles lineage information'''
    
    _delim = ';'
    _prokaryotes = ('archaea', 'bacteria')
    known_taxons = set(_prokaryotes)
    unknown = 'unknown taxonomy'
    taxonomy_re = re.compile(r'(\s|^)(\w+;((\w+\s+)*\w+;)*(\w+\s+)*\w+)(\s|$)')
    
    @classmethod
    def _register(cls, lineage):
        for t in lineage:
            cls.known_taxons.add(t)
    
    def __new__(cls, line, delimiter=';'):
        line = line.strip(' \n\r')
        taxons = [t1 for t1 in (t.strip().lower() for t in line.split(delimiter)) if t1]
        if taxons and taxons[0] in cls._prokaryotes: taxons.insert(0, 'prokaryotes')
        return tuple.__new__(cls, taxons)
    
    def __init__(self, line, delimiter=';'):
        self.str = self._delim.join(self)
        self._register(self)
    
    def __str__(self):
        return ';'.join(l.capitalize() for l in self)

    def __repr__(self): return str(self)
    
    def __eq__(self, other):
        return self.str == other.str
    
    def __sub__(self, other):
        return Lineage(self.str.replace(other.str, '', 1).strip(self._delim), 
                       delimiter=self._delim)
    
    def __getitem__(self, index):
        try: taxons = tuple.__getitem__(self, index)
        except IndexError: taxons = ''
        if isinstance(taxons, basestring): return taxons
        return self.from_iter(taxons)
    #needs to be overriden because subclassing a builtin tuple
    def __getslice__(self, i, j): return self[i:j:1]
    
    def capitalize(self, delimiter='; '):
        return delimiter.join(t.capitalize() for t in self)
    
    @property
    def first(self): return self[0]
    
    @property
    def last(self): return self[-1]
    
    def includes(self, other):
        return other.str.startswith(self.str)
    
    @classmethod
    def from_iter(cls, taxons):
        if not taxons: return cls('')
        return cls(';'.join(taxons))
    
    @classmethod
    def from_record(cls, rec):
        try: return cls.from_iter(rec.annotations['taxonomy'])
        except KeyError:
            m = Lineage.taxonomy_re.search(rec.name+' '+rec.description)
            return cls(m.group(2) if m else cls.unknown)
    
    @classmethod
    def common(cls, lineages, delimiter=';'):
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
        return all(l1 == l0 for l1 in lineages[1:])

    @classmethod
    def samefirst(cls, lineages, first=1):
        if not lineages or not first: return False
        l0 = lineages[0][:first]
        return l0 and all(l1[:first] == l0 for l1 in lineages[1:])

    @classmethod
    def samelast(cls, lineages, last=1):
        if not lineages or not last: return False
        ind = -last
        l0 = lineages[0][ind]
        return l0 and all(l1[ind] == l0 for l1 in lineages[1:])
#end class


class Organism(object):
    '''Holds an abstract organism record'''
    def __init__(self, oid, description, lineage):
        self.id = oid
        self.description = description
        self.lineage = (lineage if isinstance(lineage, Lineage) 
                        else Lineage(lineage))
        self.description = Lineage.taxonomy_re.sub('', self.description)

    def __str__(self):
        return '%s: %s [%s]' % (self.id, self.description, self.lineage)
          
    def belongs2(self, lineage):
        return lineage.includes(self.lineage)
    
    def fix_lineage(self):
        genus = self.description.split()[0].lower()
        if genus in Lineage.known_taxons:
            if genus not in self.lineage: 
                self.lineage = Lineage.from_iter(list(self.lineage)+[genus])
            else:
                self.lineage = Lineage.from_iter(self.lineage[:self.lineage.index(genus)+1])
    
    @classmethod
    def from_record(cls, rec):
        return cls(rec.id, 
                   rec.annotations.get('organism') 
                   or pretty_rec_name(rec),
                   Lineage.from_record(rec))
#end class

class Organisms(dict):
    '''Database of Organisms'''
    def add(self, rec):
        org = Organism.from_record(rec)
        if org.id in self: 
            raise KeyError('Organism "%s" is already present in database.' % org.id)
        self[org.id] = org
    
    def fix_lineages(self):
        for org in self.itervalues():
            org.fix_lineage()
            
    @property
    def common_lineage(self):
        return Lineage.common(org.lineage for org in self.itervalues())
    
    def taxons_of_rank(self, r):
        '''
        @return: a list o taxons of rank 'r' starting with 0
        '''
        taxons = []
        for org in self.itervalues():
            if len(org.lineage) > r:
                taxons.append(org.lineage[r])
        return taxons
    
    @classmethod
    def from_records(cls, records):
        orgs = cls()
        for rec in records: orgs.add(rec)
        orgs.fix_lineages()
        return orgs
    
    @classmethod
    def from_seqfiles(cls, filenames):
        db = SeqView()
        if not db.load(filenames): return None
        return cls.from_records(db)
#end class
