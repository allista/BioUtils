# coding=utf-8
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
Created on Aug 1, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import os
import tempfile
import itertools
import csv
import re
from time import time, sleep
from datetime import timedelta
from copy import deepcopy
from random import shuffle

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW, NCBIXML, Applications
from Bio import Entrez

from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from BioUtils.Tools.tmpStorage import shelf_result, roDict

from BioUtils.Tools import user_message, Progress, ProgressCounter, retry
from BioUtils.SeqUtils import mktmp_fasta, cat_records, Translator

class BatchEntrez(object):
    #defaults
    RETRIES    = 3
    PAUSE_EACH = 100
    BATCH      = 20
    PAUSE      = 60
    
    #id types
    NUC     = 'nucleotide'
    PROT    = 'protein'
    WGS     = 'wgs'
#    MGA     = 'mga'
    UPROT   = 'uniprot'
    REFSEQN = 'refseqn'
    REFSEQP = 'refseqp'
    #regexp
    ID_TYPES_RE = {
                   NUC:     re.compile(r'\b([a-zA-Z]{1}\d{5}|[a-zA-Z]{2}\d{6})\b'),
                   PROT:    re.compile(r'\b([a-zA-Z]{3}\d{5})\b'),
                   WGS:     re.compile(r'\b([a-zA-Z]{4}\d{2}\d{6,8})\b'),
#                   MGA:     re.compile(r'\b([a-zA-Z]{5}\d{7})\b'),
                   UPROT:   re.compile(r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\b'),
                   REFSEQN: re.compile(r'\b([ANX][CGTWSZMR]_\d+)\b'),
                   REFSEQP: re.compile(r'\b([ANYXZW]P_\d+)\b'),
                   }
    ID_TO_DB = {
                NUC:     'nucleotide',
                PROT:    'protein',
                WGS:     'nucleotide',
#                MGA:     '',
                UPROT:   'protein',
                REFSEQN: 'nucleotide',
                REFSEQP: 'protein',
                }
    
    def __init__(self, email):
        self.email       = email
        self._start_time = -1
    #end def
    
    @classmethod
    def extract_id(cls, id_str):
        for idt, id_re in cls.ID_TYPES_RE.items():
            match = id_re.search(id_str)
            if match is not None:
                return match.group(0), cls.ID_TO_DB[idt]
        return None, None
    #end def
    
    def _get_records(self, query, database, **efetch_kwargs):
        results = retry(lambda : Entrez.read(Entrez.esearch(db=database, term=query, usehistory="y")),
                        'Unable to get Entrez ID for query: %s' % query, self.RETRIES)
        if not results['IdList']:
            print('NCBI returned no result for query: %s' % query) 
            return []
        webenv    = results['WebEnv']
        query_key = results['QueryKey']
        #fetch genbank data for the received IDs 
        num_results = len(results['IdList'])
        print('Downloading data...')
        data = retry(lambda : Entrez.efetch(db=database, rettype="gb", retmode="text",
                                            retstart=0, retmax=num_results,
                                            webenv=webenv, query_key=query_key, **efetch_kwargs),
                     'Unable to download data for IDs: %s' % str(results['IdList']), self.RETRIES)
        #parse received data
        try: records = list(SeqIO.parse(data, 'gb'))
        except Exception as e:
            print 'Unable to parse fetched data as SeqRecords'
            print e
            return []
        finally: data.close()
        print('Done. Elapsed time: %s\n' % timedelta(seconds=time()-self._start_time))
        return records
    
    def _get_records_for_terms(self, terms, start, stop, database, **efetch_kwargs):
        query = ' or '.join(terms[i] for i in xrange(start, stop))
        num_terms = len(terms)
        print('[%3.1f%%] performing query for terms %d-%d of %d' 
              % (float(stop)/num_terms*100, start+1, stop, num_terms))
        return self._get_records(query, database, **efetch_kwargs)
    #end def
    
    def get_records(self, query, database, **efetch_kwargs):
        self._start_time = time()
        Entrez.email = self.email
        Entrez.tool = 'BatchEntrez.get_records'
        return self._get_records(query, database, **efetch_kwargs)
    
    def get_records_for_terms(self, terms, database, **efetch_kwargs):
        self._start_time = time()
        #check number of queries and warn the user
        num_terms   = len(terms)
        num_queries = num_terms/self.BATCH
        num_pauses  = 0; pause_time = 0 
        if num_queries > self.PAUSE_EACH:
            num_pauses = num_queries/self.PAUSE_EACH
            pause_time = num_pauses * self.PAUSE
            self.PAUSE_EACH = num_queries/(num_pauses+1)+1
            print('WARNING: %d separate Entrez queries will be made.\n'
                  'To comply with NCBI rules the queries will be made\n'
                  'in series of %d with %d sec pause in between.\n' 
                  % (num_queries, self.PAUSE_EACH, self.PAUSE))
            print('Total pause time will be:\n%s\n' % timedelta(seconds=pause_time))
        query_time = num_queries * 1/3.0
        if query_time > 5:
            print('No more than 3 requests per second is allowed by NCBI,\n'
                  'so *minimum* time spend for your query will be:\n%s\n' % timedelta(seconds=query_time))
        if pause_time > 0 and query_time > 5:
            print('Total *minimum* estimated time:\n%s\n' % timedelta(seconds=pause_time+query_time))
            print('Note, that depending on the load of NCBI servers it\n'
                  'may take several times as much.\n')
        #setup Entrez engine
        Entrez.email = self.email
        Entrez.tool = 'BatchEntrez.get_records_for_terms'
        #perform queries in batches
        pause_num = self.PAUSE_EACH
        records = []
        for i in xrange(0, num_terms, self.BATCH):
            if i/self.BATCH > pause_num:
                print('Pausing for %d seconds...\n' % self.PAUSE)
                sleep(self.PAUSE)
                pause_num += self.PAUSE_EACH
            try: records.extend(self._get_records_for_terms(terms, i, min(i+self.BATCH, num_terms), database, **efetch_kwargs))
            except RuntimeError as e: 
                print(e)
                continue
        print('Done.\nTotal elapsed time: %s\n' % timedelta(seconds=time()-self._start_time))
        return records
    #end def
#end class

class BlastWWW(object):
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
#            for alignment in record.alignments:
#                SQID = alignment.hit_id
#                try:
#                    entrez_ID = Entrez.esearch(db='nucleotide', term=('%s[SQID]' % SQID))
#                except Exception, e:
#                    print '\nUnable to get Entrez ID for sequence %s - %s' % (SQID, alignment.hit_def)
#                    print 'Exception message is:', e.message
#                    continue  
#                for hsp in alignment.hsps:
#                    seq_start = hsp.sbjct_start
#                    seq_end   = hsp.sbjct_end
#                    try:
#                        ehandle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', 
#                                                seq_start=str(seq_start), seq_end=str(seq_end),
#                                                id=entrez_ID)
#                        seq_record = SeqIO.read(ehandle, format='gb', alphabet=IUPAC.ambiguous_dna)
#                    except Exception, e:
#                        print '\nUnable to fetch sequence: %s - %s [%d:%d]' \
#                        % (SQID, alignment.hit_def, seq_start, seq_end)
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
    
    _fetch_targets = ('record', 'alignment', 'hsp')
    
    @classmethod
    def fetch_results(cls, email, results, from_dbs=None, what='record', **kwargs):
        if what not in cls._fetch_targets:
            raise ValueError('"what" parameter may only be in: %s' % str(cls._fetch_targets))
        if not email:
            raise ValueError('You should always provide a valid e-mail '
                             'to NCBI when performing an Entrez query.')
        #search the ID of each blast result, then fetch corresponding part of the sequence
        terms = {}
        fetched = []
        be = BatchEntrez(email)
        num_terms = len([a for r in results for a in r.alignments])
        with ProgressCounter('', num_terms, replace=False) as prg:
            for record in results:
                for alignment in record.alignments:
                    term, database = BatchEntrez.extract_id(alignment.hit_id+' '+alignment.hit_def)
                    if not term:
                        print 'Unable to extract record ID from:'
                        print alignment 
                        continue
                    if what == 'record':
                        if database in terms: 
                            terms[database].append(term)
                        else: terms[database] = [term]
                        continue
                    elif what == 'alignment':
                        loc = [0,0]
                        for hsp in alignment.hsps:
                            loc[0] = min(loc[0], hsp.sbjct_start) if loc[0] > 0 else hsp.sbjct_start
                            loc[1] = max(loc[1], hsp.sbjct_end)
                        if loc[0] > 0:
                            fetched += be.get_records(term, database,
                                                      seq_start=str(loc[0]), 
                                                      seq_end=str(loc[1]))
                    elif what == 'hsp':  
                        for hsp in alignment.hsps:
                            fetched += be.get_records(term, database,
                                                      seq_start=str(hsp.sbjct_start), 
                                                      seq_end=str(hsp.sbjct_end))
                    prg.count()
            for db in terms:
                if not from_dbs or db in from_dbs:
                    fetched += be.get_records_for_terms(terms[db], db)
        return fetched
    #end def
#end class


class BlastCLI(MultiprocessingBase):

    _clines = {
            'blastp'     : Applications.NcbiblastpCommandline, 
            # compares an amino acid query sequence against a protein sequence database
            
            'blastn'     : Applications.NcbiblastnCommandline, 
            # compares a nucleotide query sequence against a nucleotide sequence database
            
            'blastx'     : Applications.NcbiblastxCommandline, 
            # compares a nucleotide query sequence translated in all reading frames against a protein sequence database
            
            'tblastn'    : Applications.NcbitblastnCommandline, 
            # compares a protein query sequence against a nucleotide sequence database dynamically translated in all reading frames
            
            'tblastx'    : Applications.NcbitblastxCommandline, 
            # compares the six-frame translations of a nucleotide query sequence against the six-frame translations of a nucleotide sequence database. Please note that tblastx program cannot be used with the nr database on the BLAST Web page.
            
            'psiblast'   : Applications.NcbipsiblastCommandline, # Position-Specific Initiated BLAST
            'rpsblast'   : Applications.NcbirpsblastCommandline, # Reverse Position Specific BLAST
            'rpstblastn' : Applications.NcbirpstblastnCommandline, # Translated Reverse Position Specific BLAST
            'deltablast' : Applications.NcbideltablastCommandline # Protein-Protein domain enhanced lookup time accelerated blast
            }

    def __init__(self, abort_event):
        super(BlastCLI, self).__init__(abort_event)
    
    @classmethod
    def blast(cls, command, **kwargs):
        '''Generic wrapper for blast commandline programs'''
        results_file = kwargs.pop('save_results_to', None)
        parse_results = kwargs.pop('parse_results', False)
        if results_file: bout = results_file
        else:
            parse_results = True
            f, bout = tempfile.mkstemp('.xml')
            os.close(f)
        try:
            cline = cls._clines[command]
            cmd = cline(outfmt=5, out=bout, **kwargs)
            out, err = cmd()
            if err:
                print '\nError while performing local %s:' % command
                print cmd
                print out
                print err
                return None
            if parse_results:
                with open(bout) as inp:
                    results = list(NCBIXML.parse(inp))
                return results
        except Exception, e:
            print '\nError while performing local %s:' % command
            print e
            return None
        finally:
            if not results_file: os.unlink(bout)
        return None
    
    @classmethod
    def blast_seq(cls, query, db='nr', evalue=0.001, command='blastn', **kwargs):
        '''Perform local blast of a SeqRecord against a database'''
        qfile = mktmp_fasta(query)
        results = cls.blast(command, query=qfile, db=db, 
                                evalue=evalue, **kwargs)
        os.unlink(qfile)
        return results

    @classmethod
    def s2f_blast(cls, query, subject_file, evalue=0.001, command='blastn', **kwargs):
        '''Perform local blast of a SeqRecord against a file with sequence(s)'''
        qfile = mktmp_fasta(query)
        results = cls.blast(command, query=qfile, subject=subject_file, 
                            evalue=evalue, **kwargs)
        os.unlink(qfile)
        return results
    
    @classmethod
    def s2s_blast(cls, query, subject, evalue=0.001, command='blastn', **kwargs):
        '''Perform local blast of one SeqRecord against the other'''
        sfile = mktmp_fasta(subject)
        results = cls.s2f_blast(query, sfile, evalue, command, **kwargs)
        os.unlink(sfile)
        return results
    
    def s2s_blast_batch(self, queries, subjects, subject_locs=None, evalue=0.001, command='blastn', **kwargs):
        results = self._s2s_blast_batch(queries, subjects, subject_locs, evalue, command, **kwargs)
        if not results: return None
        with ProgressCounter('Parsing results...', len(results)) as prg:
            for qi in xrange(len(results)):
                for si in xrange(len(results[qi])):
                    results_name = results[qi][si]
                    if not results_name: continue
                    with roDict(results_name) as db:
                        results[qi][si] = db['result']
                prg.count()
        return results
    
    @classmethod
    def all_hsps(cls, blast_results, max_rlen=0):
        if not blast_results: return None
        hsps = []
        for r in blast_results:
            for alignment in r.alignments:
                if max_rlen > 0:
                    hsps += [hsp for hsp in alignment.hsps 
                             if float(len(hsp.query))/r.query_length >= max_rlen]
                else: hsps += alignment.hsps
        return hsps or None
    
    def _s2s_blast_batch(self, queries, subjects, subject_locs=None, evalue=0.001, command='blastn', **kwargs):
        queries_len = len(queries)
        subjects_len = len(subjects)
        results = [[None for _s in subjects] for _q in queries]
        pairs = list(itertools.product(xrange(queries_len), xrange(subjects_len)))
        ignore_none_locs = kwargs.pop('ignore_none_locs', False)
        shuffle(pairs)
        @MultiprocessingBase.data_mapper
        @shelf_result
        def _worker(qs, queries, subjects, subject_locs):
            query = queries[qs[0]]
            subject = subjects[qs[1]]
            if query is None or subject is None: return None
            if subject_locs:
                if subject_locs[qs[1]]: 
                    loc = tuple(subject_locs[qs[1]])
                    kwargs['subject_loc'] = '%d-%d' % loc
                elif ignore_none_locs: return None
                elif 'subject_loc' in kwargs: 
                    del kwargs['subject_loc']
            elif 'subject_loc' in kwargs: del kwargs['subject_loc']
            return BlastCLI.s2s_blast(query, subject, evalue, command, **kwargs)
        @MultiprocessingBase.results_assembler
        def _assembler(index, hsps, results, pairs, prg):
            qs = pairs[index]
            results[qs[0]][qs[1]] = hsps
            prg.count()
        with ProgressCounter('Performing multiple %s searches:'%command, len(pairs)) as prg:
            work = self.Work()
            work.prepare_jobs(_worker, pairs, None, queries, subjects, subject_locs)
            work.set_assembler(_assembler, results, pairs, prg)
            self.start_work(work)
            if not self.wait(work): return None
        return results
    
    @staticmethod
    @MultiprocessingBase.data_mapper
    def _find_features_by_hsps(qs, translations, blast_results):
        results_name = blast_results[qs[0]][qs[1]]
        if not results_name: return None
        with roDict(results_name) as db:
            hsps = BlastCLI.all_hsps(db['result'])
        if not hsps:
            print 'no blast results in %s: %s' % (results_name, hsps)  
            return None
        hsp = min(hsps, key=lambda h: h.expect)
        for f in translations[qs[1]].features:
            if hsp.sbjct_start-1 in f and hsp.sbjct_end-1 in f:
                aln_len = abs(hsp.sbjct_start-hsp.sbjct_end)+1
                return f, aln_len, float(hsp.identities)/hsp.align_length, hsp.expect
        return None
        
    @staticmethod
    def _features_of_type_i(record, ftype):
        return [i for i, f in enumerate(record.features) if f.type == ftype]
        
    @staticmethod
    def _get_genes(rec):
        features = BlastCLI._features_of_type_i(rec, 'CDS')
        if not features:
            features = BlastCLI._features_of_type_i(rec, 'gene')
            if not features:
                print 'No gene/CDS features found in:'
                print rec.id, rec.description
                return None
        return features
    
    def _get_fois(self, records, foi):
        for q in foi: foi[q] = re.compile(foi[q])
        num_records = len(records)
        def _get_foi(ri, records, foi):
            r = records[ri]
            features = []
            for fi, f in enumerate(r.features):
                for q in foi:
                    qv = f.qualifiers.get(q)
                    if qv and foi[q].search(qv[0]):
                        features.append(fi)
            for fi, f in enumerate(features):
                srec = r.features[f].extract(r)
                genes = self._get_genes(srec)
                if genes: features[fi] = [srec.features[gi].qualifiers.get('gene_id') for gi in genes]
            return features
        return self.parallelize_work(1, _get_foi, range(num_records), records, foi)
    
    def g2g_blastp(self, reference, subjects, table='Standard', 
                   evalue=0.001, max_rlen=0, features_of_interest=None):
        '''
        Perform blastp of each coding sequence of the reference against each 
        subject, which is first translated gene-by-gene.
        Parameters
        @param reference: SeqRecord object of the reference genome
        @param subjects: a list of SeqRecord objects of subject genomes
        @param table: translation table number (see NCBI site for description)
        @param evalue: filter out blastp results with E-value grater than this
        @param max_rlen: filter out blastp results which are shorter than this 
        fraction of target gene length
        @param features_of_interest: list of dictionaries of the form 
        {qualifier_name : qualifier_value}
        to mark features denoting known clusters that should be analyzed one 
        against the other
        @return: list of pairs (CDS, (blast_result1, blast_result2, ...)) 
        where CDS is a gene/CDS feature from the reference.features list 
        and blast_resultN is a list of results for the N-th  
        subject, containing following information:
        (hit_feature, align_length, percent_identity, evalue)
        where hit_feature is a SeqFeature object of the gene/CDS of the subject
        where top blast hit is located, align_length is the length of the hit,
        percent_identity is the ratio of number of identities and align_length [0; 1]
        and evalue is the E-value of the top hit.
        '''
        if not reference or not subjects:
            print 'No reference or subject sequences provided' 
            return None
        #get list of features to query
        with user_message('Searching for gene/CDS features in provided sequences...'):
            all_records = [reference]+subjects
            num_records = len(all_records)
            features = self.parallelize_work(1, lambda ri, records: self._get_genes(records[ri]), 
                                             range(num_records), 
                                             all_records)
            if self.aborted():
                print '\nAborted'
                return None
            if not features or not features[0]:
                print '\nReference sequence does not contain annotated genes:'
                print reference.id, reference.description
                return None
            if len([f for f in features if f]) < 2:
                print '\nSubject sequences do not contain annotated genes'
                return None
            #add gene ids
            for ri, genes in enumerate(features):
                if not genes: continue
                r = all_records[ri]
                for gene_id, gi in enumerate(genes):
                    r.features[gi].qualifiers['feature_id'] = gi
                    r.features[gi].qualifiers['gene_id'] = gene_id
        #get features of interest if requested
        fois = None
        if features_of_interest:
            with user_message('Searching for features of interest...'):
                fois = []
                for foi in features_of_interest:
                    foi = self._get_fois(all_records, foi)
                    if foi and foi[0]: fois.append(foi)
                    if self.aborted():
                        print '\nAborted'
                        return None
        #translate features to proteins
        with Progress('Translating genes found in the reference and subjects...', num_records) as prg:
            translator = Translator(self._abort_event)
            translations = [None]*num_records
            foi_translations = [[None]*num_records for _f in fois]
            for i, (f, rec) in enumerate(zip(features, all_records)):
                if not f:
                    prg.step(i) 
                    continue
                translation = translator.translate(rec, f, table)
                if not translation: return None 
                if i > 0: 
                    translations[i] = cat_records(translation)
                    if fois:
                        for ifoi, foi in enumerate(fois):
                            foi_loc = [0, 0]
                            for foi_var in foi[i]: 
                                if not foi_var: continue
                                for gid in foi_var:
                                    l = translations[i].features[gid].location
                                    foi_loc[0] = min(int(l.start)+1, foi_loc[0]) if foi_loc[0] > 0 else int(l.start)+1
                                    foi_loc[1] = max(int(l.end), foi_loc[1])
                            if foi_loc[0] > 0: foi_translations[ifoi][i] = foi_loc 
                else: 
                    translations[i] = translation
                    if fois: 
                        for ifoi, foi in enumerate(fois):
                            foi_translations[ifoi][i] = [[translation[gid] for gid in foi_var] for foi_var in foi[i]]
                prg.step(i)
        #blast features against subjects
        with user_message('Performing local blast of every translated gene in the reference against every translated subject...', '\n'):
            stranslations = translations[1:]
            blast_results = self._s2s_blast_batch(translations[0], stranslations, None, evalue, 
                                                  command='blastp', task='blastp')
            if self.aborted():
                print '\nAborted'
                return None
            if not blast_results:
                print '\nBlast have not returned any results.' 
                return None
        if fois: #redo blast for fois and replace the results
            with user_message('Rerunning blast for FOIs...', '\n'):
                for ifoi, foi in enumerate(foi_translations):
                    sfoi_locs = foi[1:]
                    for i, foi_var in enumerate(foi[0]):
                        foi_blast = self._s2s_blast_batch(foi_var, stranslations, sfoi_locs, evalue, 
                                                          command='blastp', task='blastp')
                        if self.aborted():
                            print '\nAborted'
                            return None
                        if not foi_blast: continue
                        for gi, gid in enumerate(fois[ifoi][0][i]):
                            if foi_blast[gi]:
                                blast_results[gid] = foi_blast[gi]
        #process blast results
        pairs = list(itertools.product(xrange(len(translations[0])), xrange(len(stranslations))))
        with ProgressCounter('Searching for genes in subjects that overlap with top blast hits...', len(pairs)) as prg:
            work = self.Work()
            work.prepare_jobs(self._find_features_by_hsps, pairs,
                              None, stranslations, blast_results)
            @MultiprocessingBase.results_assembler
            def _assembler(index, result, blast_results, pairs, prg):
                qs = pairs[index]
                blast_results[qs[0]][qs[1]] = result
                prg.count()
            work.set_assembler(_assembler, blast_results, pairs, prg)
            self.start_work(work)
            if not self.wait(work): return None
        return zip((reference.features[f] for f in features[0]), blast_results)
    
    @staticmethod
    def g2g_to_csv(filename, reference, subjects, g2g_results):
        locus_tag = lambda g, f, i: f.qualifiers.get('locus_tag', ['%s_%s' % (g.name, i+1)])[0]
        product = lambda f: f.qualifiers.get('product', [''])[0]
        with open(filename, 'wb') as out:
            writer = csv.writer(out)
            header = [reference.name+'_gene', 'locus_tag', reference.name+'_desc']
            for g in subjects:
                header += [g.name+'_gene', 'locus_tag', g.name+'_desc', g.name+'_coverage', g.name+'_percent', g.name+'_evalue']
            writer.writerow(header)
            for i, (feature, res) in enumerate(g2g_results):
                if not res: continue
                row = [i+1, locus_tag(reference, feature, i), product(feature)]
                for sres, g in zip(res, subjects):
                    if sres is None:
                        row += ['', '', '', '', '', '']
                    else:
                        f, align_length, percent_identity, evalue = sres
                        gid = f.qualifiers.get('gene_id', 'N/A')
                        row += [gid, locus_tag(g, f, gid), product(f), 
                                '%d/%d' % (align_length, len(f)), percent_identity*100, evalue]
                writer.writerow(row)
#end class