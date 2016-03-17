# coding=utf-8

'''
Created on Mar 17, 2016

@author: Allis Tauri <allista@gmail.com>
'''

from BioUtils.Tools.Output import ProgressCounter
from BioUtils.Tools.Misc import ListDB

from .BatchEntrez import BatchEntrez
from .BlastBase import _BlastBase

class BlastWWW(_BlastBase):
    '''Perform a blast search possibly restricted by Entrez query, then fetch
    results as sequences with annotations'''
    
#    def __init__(self, job_id):
#        self._job_id = job_id
#        self._sequences_filename = self._job_id + '-blast-seqs.gb'
#        self._results_filename   = self._job_id + '-blast.xml'
#        self._blast_results = None
#        self._blast_results_sequences = None
#        self.have_results   = False
#    #end def
#    
#    def blast(self, query, command='blastn', expect=10, megablast=False, entrez_query='', **kwargs):
#        if not query: return False
#        print '\nStarting a BLAST search. This may take awhile...'
#        try:
#            blast_results = NCBIWWW.qblast(command, 
#                                           self.database, 
#                                           self._query.format('fasta'),
#                                           expect       = expect, 
#                                           megablast    = megablast, 
#                                           entrez_query = entrez_query,
#                                           **kwargs)
#            #save results to a file
#            results_file = open(self._results_filename, 'w')
#            results_file.write(blast_results.read())
#            results_file.close()
#            blast_results.close()
#            print '\nBlast output was written to:\n   ' + self._results_filename
#            #parse results
#            results_file  = open(self._results_filename, 'r')
#            self._blast_results = list(NCBIXML.parse(results_file))
#            results_file.close()
#            #clear fetched sequences
#            self._blast_results_sequences = None
#        except Exception, e:
#            print 'BlastWWW.blast: failed to obtain BLAST search results from NCBI.'
#            print e
#            return False
#        self.have_results = True
#    #end def
#    
#    
#    def load_results(self):
#        #load blast results
#        print '\nLoading previously saved BLAST results:\n   ', self._results_filename
#        try:
#            results_file        = open(self._results_filename, 'r')
#            self._blast_results = list(NCBIXML.parse(results_file))
#            results_file.close()
#            self.have_results = True
#        except Exception, e:
#            print '\nFailed to load blast results.\n', e
#    #end def
#    
#    
#    def fetch_results_(self, email):
#        if not self.have_results: return
#        #if sequences were already fetched, just return them
#        if self._blast_results_sequences:
#            return deepcopy(self._blast_results_sequences)
#        if not email:
#            raise ValueError('\nYou should always provide a valid e-mail '
#                             'to NCBI when performing an Entrez query.')
#        print '\nFetching records from Entrez database. This may take awhile...'
#        Entrez.email = email
#        #search the ID of each blast result, then fetch corresponding part of the sequence
#        self._blast_results_sequences = list()
#        for record in self._blast_results:
#            for region in record.alignments:
#                SQID = region.hit_id
#                try:
#                    entrez_ID = Entrez.esearch(db='nucleotide', term=('%s[SQID]' % SQID))
#                except Exception, e:
#                    print '\nUnable to get Entrez ID for sequence %s - %s' % (SQID, region.hit_def)
#                    print 'Exception message is:', e.message
#                    continue  
#                for hsp in region.hsps:
#                    seq_start = hsp.sbjct_start
#                    seq_end   = hsp.sbjct_end
#                    try:
#                        ehandle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', 
#                                                seq_start=str(seq_start), seq_end=str(seq_end),
#                                                id=entrez_ID)
#                        seq_record = SeqIO.read(ehandle, format='gb', alphabet=IUPAC.ambiguous_dna)
#                    except Exception, e:
#                        print '\nUnable to fetch sequence: %s - %s [%d:%d]' \
#                        % (SQID, region.hit_def, seq_start, seq_end)
#                        print e
#                        continue
#                    self._blast_results_sequences.append(seq_record)
#        #write fetched sequences to a file
#        if self._blast_results_sequences:
#            try:
#                sequences_file      = open(self._sequences_filename, 'w')
#                SeqIO.write(self._blast_results_sequences, sequences_file, 'gb')
#                sequences_file.close()
#            except Exception, e:
#                print '\nFailed to save fetched sequences.'
#                print  e
#        return deepcopy(self._blast_results_sequences)
#    #end def
#
#    def load_sequences(self):
#        #load fetched sequences
#        print '\nLoading previously saved sequences:\n   ', self._sequences_filename
#        try:
#            sequences_file      = open(self._sequences_filename, 'r')
#            self._blast_results = list(SeqIO.parse(sequences_file, 'gb', IUPAC.ambiguous_dna))
#            sequences_file.close()
#        except Exception:
#            print '\nFailed to load fetched sequences.'
#    #end def

    @classmethod
    def _fetch_query(cls, q, entrez):
        if q.region:
            records = entrez.get_records(q.term, q.db,
                                         seq_start=q.region.start,
                                         seq_end=q.region.end,
                                         strand=q.entrez_strand)
            return q.expand_subregions(records)
        else: return entrez.get_records(q.term, q.db)
    
    @classmethod
    def fetch_results(cls, email, results, from_dbs=None, what='record', **kwargs):
        if what not in cls._fetch_targets:
            raise ValueError('"what" parameter may only be in: %s' % str(cls._fetch_targets))
        if not email:
            raise ValueError('You should always provide a valid e-mail '
                             'to NCBI when performing an Entrez query.')
        #search the ID of each blast result, then fetch corresponding part of the sequence
        batch = ListDB()
        single = ListDB()
        fetched = []
        be = BatchEntrez(email)
        with ProgressCounter('', 0, replace=False) as prg:
            for record in results:
                for alignment in record.alignments:
                    q = cls.Query(alignment, what)
                    if not q: continue
                    if q.region: single[q.db] = q
                    else: batch[q.db] = q
                    prg.count()
            for db in single:
                if not from_dbs or db in from_dbs:
                    for q in single[db]: 
                        fetched += cls._fetch_query(q, be)
                        prg.count()
            for db in batch:
                if not from_dbs or db in from_dbs:
                    fetched += be.get_records_for_terms(batch[db], db)
                    prg.cound()
        return fetched
    #end def
#end class
