# coding=utf-8

'''
Created on Mar 17, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import re
import inspect
from collections import Counter, OrderedDict

def num_alignments(results):
    return sum(len(rec.alignments) for rec in results)

def print_hsps(results):
    for rec in results:
        for ali in rec.alignments:
            for hsp in ali.hsps:
                print hsp
                print


class _BlastBase(object):
    '''Base class for both standalone and WWW versions'''
    
    _fetch_targets = ('record', 'alignment', 'hsp')
    
    class Query(object):
        class Region(object):
            def __init__(self, start, end):
                se = sorted((start, end))
                self.start = se[0]
                self.end = se[1]
                self.strand = start < end
                
            def subrec(self, rec):
                srec = rec[self.start:self.end]
                print '='*80
                return srec if self.strand else srec.reverse_complement()
            
            @property
            def entrez_strand(self):
                return 1 if self.strand else 2
            
            @property
            def blast_strand(self):
                return 'plus' if self.strand else 'minus'
            
            def __str__(self):
                return '%d-%d %s' % (self.start, self.end, self.blast_strand)
        
        def __init__(self, alignment, what='record'):
            self.region = None
            self.subregions = []
            self.term, self.db = BlastID.extract(alignment.hit_id+' '+alignment.hit_def)
            if what == 'record': return
            if len(alignment.hsps) == 1:
                hsp = alignment.hsps[0] 
                self.region = self.Region(hsp.sbjct_start, hsp.sbjct_end)
                return
            loc = [0,0]
            for hsp in alignment.hsps:
                self.subregions.append(self.Region(hsp.sbjct_start, hsp.sbjct_end))
                loc[0] = (min(loc[0], hsp.sbjct_start, hsp.sbjct_end) if loc[0] > 0 
                          else min(hsp.sbjct_start, hsp.sbjct_end))
                loc[1] = max(loc[1], hsp.sbjct_start, hsp.sbjct_end)
            if what != 'hsp':
                strand = Counter(r.strand for r in self.subregions)
                if strand.most_common(1)[0]: self.region = self.Region(*loc)
                else: self.region = self.Region(loc[1], loc[0])
                self.subregions = None
            else: self.region = self.Region(*loc)
            
        def __nonzero__(self): return bool(self.term)
        
        def expand_subregions(self, records):
            if not self.subregions: return records
            subrecords = []
            for rec in records:
                subrecords.extend(r.subrec(rec) for r in self.subregions)
            return subrecords
        
        def __str__(self):
            if self.subregions:
                return '\n'.join('%s %s' % (self.term, r) for r in self.subregions)
            if self.region: 
                return '%s %s' % (self.term, self.region)
            return self.term
#end class

class BlastID(object):
    #id types
    NUC     = 'nucleotide'
    PROT    = 'protein'
    WGS     = 'wgs'
#    MGA     = 'mga'
    UPROT   = 'uniprot'
    REFSEQN = 'refseqn'
    REFSEQP = 'refseqp'
    LOCAL   = '0local'
    #regexp
    ID_TYPES_RE = OrderedDict(((LOCAL,   re.compile(r'\b(gnl\|BL_ORD_ID\|\d+)\b')),
                               (NUC,     re.compile(r'\b([a-zA-Z]{1}\d{5}(\.\d+)*|[a-zA-Z]{2}\d{6}(\.\d+)*)\b')),
                               (PROT,    re.compile(r'\b([a-zA-Z]{3}\d{5}(\.\d+)*)\b')),
                               (WGS,     re.compile(r'\b([a-zA-Z]{4}\d{2}\d{6,8}(\.\d+)*)\b')),
                               #(MGA,     re.compile(r'\b([a-zA-Z]{5}\d{7})\b')),
                               (UPROT,   re.compile(r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(\.\d+)*)\b')),
                               (REFSEQN, re.compile(r'\b([ANX][CGTWSZMR]_\d+(\.\d+)*)\b')),
                               (REFSEQP, re.compile(r'\b([ANYXZW]P_\d+(\.\d+)*)\b'))))

    ID_TO_DB = {
                NUC:     'nucleotide',
                PROT:    'protein',
                WGS:     'nucleotide',
#                MGA:     'mga',
                UPROT:   'protein',
                REFSEQN: 'nucleotide',
                REFSEQP: 'protein',
                LOCAL:   'local',
                }
    
    @classmethod
    def extract(cls, id_str):
        for idt, id_re in cls.ID_TYPES_RE.items():
            match = id_re.search(id_str)
            if match is not None:
                return match.group(0), cls.ID_TO_DB[idt]
        return None, None
    #end def
#end class


class BlastFilter(object):
    def __init__(self, predicate):
        self.P   = predicate
        self.AND = None
        self.OR  = None
    
    def test(self, alignment):
        return (self.P(alignment) 
                and (self.AND.test(alignment) if self.AND else True)
                or  (self.OR.test(alignment) if self.OR else False))
    
    def __call__(self, results):
        for record in results:
            for i in xrange(len(record.alignments)-1,-1,-1):
                if not self.test(record.alignments[i]):
                    del record.alignments[i]
                    
    def __str__(self):
        try: s = inspect.getsource(self.P).strip()
        except: s = 'Unknown Predicate'
        if self.AND:
            s += '\n'
            s += '\n'.join('    AND %s' % l for l in str(self.AND).splitlines())
        if self.OR: 
            s += '\nOR %s' % self.OR
        return s
