# coding=utf-8

'''
Created on Mar 17, 2016

@author: Allis Tauri <allista@gmail.com>
'''

from Bio.Blast import NCBIXML

from BioUtils.Tools.Output import ProgressCounter
from BioUtils.Tools.Misc import ListDB

from .BatchEntrez import BatchEntrez
from .BlastBase import _BlastBase

class BlastWWW(_BlastBase):
    '''Perform a blast search possibly restricted by Entrez query, then fetch
    results as sequences with annotations'''

    @classmethod
    def load_results(cls, filename):
        try:
            with open(filename) as inp:
                results = list(NCBIXML.parse(inp))
            return results
        except Exception as e:
            print 'Unhandled exception: %s' % str(e)
            return None

    @classmethod
    def _fetch_query(cls, q, entrez):
        if q.subregions:
            records = []
            for r in q.subregions:
                print 'Fetching: %s %s' % (q.term, r)
                records.extend(entrez.get_records(q.term, q.db,
                                                  seq_start=r.start,
                                                  seq_stop=r.end,
                                                  strand=r.entrez_strand))
            return records
        elif q.region:
            print 'Fetching: %s %s' % (q.term, q.region)
            return entrez.get_records(q.term, q.db,
                                      seq_start=q.region.start,
                                      seq_stop=q.region.end,
                                      strand=q.region.entrez_strand)
        else:
            print 'Fetching: %s' % q.term
            return entrez.get_records(q.term, q.db)

    def fetch_queries(self, email, queries, from_dbs=None, **kwargs):
        if not email:
            raise ValueError('You should always provide a valid e-mail '
                             'to NCBI when performing an Entrez query.')
        # search the ID of each blast result, then fetch corresponding part of the sequence
        batch = ListDB()
        single = ListDB()
        fetched = []
        entrez = BatchEntrez(self._abort_event, email)
        with ProgressCounter('', 0, replace=False) as prg:
            for q in queries:
                if self.aborted(): return []
                if q.region: single[q.db] = q
                else: batch[q.db] = q
            for db in single:
                if not from_dbs or db in from_dbs:
                    for q in single[db]:
                        if self.aborted(): return []
                        fetched += self._fetch_query(q, entrez)
                        prg.count()
            for db in batch:
                if self.aborted(): return []
                if not from_dbs or db in from_dbs:
                    fetched += entrez.get_records_for_terms([q.term for q in batch[db]], db)
                    prg.count()
        return fetched

    def fetch_results(self, email, results, from_dbs=None, what='record', **kwargs):
        if what not in self._fetch_targets:
            raise ValueError('"what" parameter may only be one of: %s' % str(self._fetch_targets))
        queries = []
        for ali in self.iter_alignments(results):
            q = self.Query(ali, what)
            if q: queries.append(q)
        return self.fetch_queries(email, queries, from_dbs, **kwargs)
#end class
