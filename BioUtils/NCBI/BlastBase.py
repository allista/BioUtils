# coding=utf-8

'''
Created on Mar 17, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import re
import inspect
from collections import Counter, OrderedDict


class _BlastBase(object):
    '''Base class for both standalone and WWW versions'''
    
    _fetch_targets = ('record', 'alignment', 'hsp')

    @staticmethod
    def num_alignments(results):
        return sum(len(rec.alignments) for rec in results)

    @staticmethod
    def have_alignments(results):
        return any(len(rec.alignments) > 0 for rec in results)

    @classmethod
    def print_hsps(cls, results):
        for hsp in cls.iter_hsps(results):
            print '%s\n' % str(hsp)

    @staticmethod
    def iter_alignments(results):
        for rec in results:
            for ali in rec.alignments: yield ali

    @classmethod
    def iter_hsps(cls, results):
        for ali in cls.iter_alignments(results):
            for hsp in ali.hsps: yield hsp

    @classmethod
    def all_hsps(cls, blast_results, min_rlen=0):
        if not blast_results: return None
        hsps = []
        for r in blast_results:
            for alignment in r.alignments:
                if min_rlen > 0:
                    hsps += [hsp for hsp in alignment.hsps
                             if float(len(hsp.query)) / r.query_length >= min_rlen]
                else: hsps += alignment.hsps
        return hsps or None

    class Query(object):
        class Region(object):
            def __init__(self, start, end):
                se = sorted((start, end))
                self.start = se[0]
                self.end = se[1]
                self.strand = start < end
                
            def subrec(self, rec):
                srec = rec[self.start:self.end]
                return srec if self.strand else srec.reverse_complement()
            
            @property
            def entrez_strand(self):
                return 1 if self.strand else 2
            
            @property
            def blast_strand(self):
                return 'plus' if self.strand else 'minus'
            
            def __str__(self):
                return '%d-%d %s' % (self.start, self.end, self.blast_strand)
        
        def __init__(self, alignment, what='record', start_offset=0, end_offset=0):
            self.region = None
            self.subregions = []
            self.term, self.db = BlastID.extract(alignment.hit_id+' '+alignment.hit_def)
            if what == 'record': return
            loc = [0,0]
            for hsp in alignment.hsps:
                if hsp.sbjct_start < hsp.sbjct_end:
                    s = max(hsp.sbjct_start-start_offset, 1)
                    e = hsp.sbjct_end+end_offset
                else:
                    s = hsp.sbjct_start + start_offset
                    e = max(hsp.sbjct_end - end_offset, 1)
                self.subregions.append(self.Region(s, e))
                loc[0] = min(loc[0], s, e) if loc[0] > 0 else min(s, e)
                loc[1] = max(loc[1], s, e)
            if what != 'hsp':
                strand = Counter(r.strand for r in self.subregions)
                if not strand.most_common(1)[0]: loc = [loc[1], loc[0]]
                self.subregions = None
            self.region = self.Region(*loc)
            
        def __nonzero__(self): return bool(self.term)
        
        def expand_subregions(self, records):
            if not self.subregions: return records
            subrecords = []
            for rec in records:
                subrecords.extend([r.subrec(rec) for r in self.subregions])
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
    def __init__(self, predicate, filter_hsps = False):
        self.filter_hsps = filter_hsps
        self.P   = predicate
        self.AND = None
        self.OR  = None
    
    def test(self, obj):
        return (self.P(obj)
                and (self.AND.test(obj) if self.AND else True)
                or (self.OR.test(obj) if self.OR else False))
    
    def __call__(self, results):
        for record in results:
            for i in xrange(len(record.alignments)-1,-1,-1):
                if self.filter_hsps:
                    for j in xrange(len(record.alignments[i].hsps) - 1, -1, -1):
                        if not self.test(record.alignments[i].hsps[j]):
                            del record.alignments[i].hsps[j]
                    if not record.alignments[i].hsps:
                        del record.alignments[i]
                elif not self.test(record.alignments[i]):
                    del record.alignments[i]
        return results
                    
    def __str__(self):
        try: s = inspect.getsource(self.P).strip()
        except: s = 'Unknown Predicate'
        if self.AND:
            s += '\n'
            s += '\n'.join('    AND %s' % l for l in str(self.AND).splitlines())
        if self.OR: 
            s += '\nOR %s' % self.OR
        return s

    def AndFilter(self, predicate, filter_hsps = False):
        if self.AND is None:
            self.AND = BlastFilter(predicate, filter_hsps)
        else: self.AND.AndFilter(predicate, filter_hsps)

    def OrFilter(self, predicate, filter_hsps = False):
        if self.OR is None:
            self.OR = BlastFilter(predicate, filter_hsps)
        else: self.OR.AndFilter(predicate, filter_hsps)


class UniqueIDs(BlastFilter):
    def __init__(self, ids=[]):
        self._ids = set(ids)
        BlastFilter.__init__(self, self._check_id, False)

    def _check_id(self, ali):
        if ali.hit_id in self._ids: return False
        self._ids.add(ali.hit_id)
        return True
