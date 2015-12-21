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
from copy import deepcopy

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW, NCBIXML, Applications
from Bio import Entrez

from DegenPrimer.MultiprocessingBase import MultiprocessingBase
from DegenPrimer.tmpStorage import shelf_result, roDict, cleanup_file

from BioUtils.CommonTools import user_message

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

class BlastFetcher(object):
    '''Make a blast search possibly restricted by Entrez query, then fetch
    results as sequences with annotations'''

    def __init__(self, job_id):
        self._job_id = job_id
        self._sequences_filename = self._job_id + '-blast-seqs.gb'
        self._results_filename   = self._job_id + '-blast.xml'
        self._blast_results = None
        self._blast_results_sequences = None
        self.have_results   = False
    #end def
    
    def blast(self, query, command='blastn', expect=10, megablast=False, entrez_query=''):
        if not query: return False
        print '\nStarting a BLAST search. This may take awhile...'
        try:
            blast_results = NCBIWWW.qblast(command, 
                                           self.database, 
                                           self._query.format('fasta'),
                                           expect       = expect, 
                                           megablast    = megablast, 
                                           entrez_query = entrez_query)
            #save results to a file
            results_file = open(self._results_filename, 'w')
            results_file.write(blast_results.read())
            results_file.close()
            blast_results.close()
            print '\nBlast output was written to:\n   ' + self._results_filename
            #parse results
            results_file  = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            #clear fetched sequences
            self._blast_results_sequences = None
        except Exception, e:
            print 'BlastFetcher.blast: failed to obtain BLAST search results from NCBI.'
            print e
            return False
        self.have_results = True
    #end def
    
    
    def load_results(self):
        #load blast results
        print '\nLoading previously saved BLAST results:\n   ', self._results_filename
        try:
            results_file        = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            self.have_results = True
        except Exception, e:
            print '\nFailed to load blast results.\n', e
    #end def
    
    
    def fetch_results(self, email):
        if not self.have_results: return
        #if sequences were already fetched, just return them
        if self._blast_results_sequences:
            return deepcopy(self._blast_results_sequences)
        if not email:
            raise Exception('\nYou should always provide a valid e-mail '
                            'to NCBI when performing an Entrez search.')
        print '\nFetching records from Entrez database. This may take awhile...'
        Entrez.email = email
        #search the ID of each blast result, then fetch corresponding part of the sequence
        self._blast_results_sequences = list()
        for record in self._blast_results:
            for alignment in record.alignments:
                SQID = alignment.hit_id
                try:
                    entrez_ID = Entrez.esearch(db='nucleotide', term=('%s[SQID]' % SQID))
                except Exception, e:
                    print '\nUnable to get Entrez ID for sequence %s - %s' % (SQID, alignment.hit_def)
                    print 'Exception message is:', e.message
                    continue  
                for hsp in alignment.hsps:
                    seq_start = hsp.sbjct_start
                    seq_end   = hsp.sbjct_end
                    try:
                        ehandle   = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', 
                                                  seq_start=str(seq_start), seq_end=str(seq_end),
                                                  id=entrez_ID)
                        seq_record = SeqIO.read(ehandle, format='gb', alphabet=IUPAC.ambiguous_dna).next()
                    except Exception, e:
                        print '\nUnable to fetch sequence: %s - %s [%d:%d]' \
                        % (SQID, alignment.hit_def, seq_start, seq_end)
                        print e
                        continue
                    self._blast_results_sequences.append(seq_record)
        #write fetched sequences to a file
        if self._blast_results_sequences:
            try:
                sequences_file      = open(self._sequences_filename, 'w')
                SeqIO.write(self._blast_results_sequences, sequences_file, 'gb')
                sequences_file.close()
            except Exception, e:
                print '\nFailed to save fetched sequences.'
                print  e
        return deepcopy(self._blast_results_sequences)
    #end def
    
    
    def load_sequences(self):
        #load fetched sequences
        print '\nLoading previously saved sequences:\n   ', self._sequences_filename
        try:
            sequences_file      = open(self._sequences_filename, 'r')
            self._blast_results = list(SeqIO.parse(sequences_file, 'gb', IUPAC.ambiguous_dna))
            sequences_file.close()
        except Exception:
            print '\nFailed to load fetched sequences.'
    #end def
#end class


def mktmp_fasta(rec):
    fd, fn = tempfile.mkstemp('.fas', 'wb')
    f = os.fdopen(fd, 'wb')
    SeqIO.write(rec, f, 'fasta')
    f.close()
    return fn

class LocalBlast(MultiprocessingBase):

    def __init__(self, abort_event):
        super(LocalBlast, self).__init__(abort_event)

    @staticmethod
    def r2r_blast(query, subject, evalue=0.001, max_rlen=0, command='blastn', **kwargs):
        '''Perform local blast of one SeqRecord against the other'''
        qfile = mktmp_fasta(query)
        sfile = mktmp_fasta(subject)
        f, bout = tempfile.mkstemp('.xml')
        os.close(f)
        try:
            cline = _clines[command]
            cmd = cline(query=qfile, subject=sfile, 
                        evalue=evalue, outfmt=5, out=bout,
                        **kwargs)
            out, err = cmd()
            if err:
                print 'Error while performing local blast:'
                print cmd
                print out
                print err
                return None
            with open(bout) as inp:
                results = list(NCBIXML.parse(inp))
            hsps = []
            qlen = len(query)
            for r in results:
                for alignment in r.alignments:
                    if max_rlen > 0:
                        hsps += [hsp for hsp in alignment.hsps 
                                 if float(len(hsp.query))/qlen >= max_rlen]
                    else: hsps += alignment.hsps
            return hsps or None
        finally:
            os.unlink(qfile)
            os.unlink(sfile)
            os.unlink(bout)
        return None
    
    def r2r_blast_batch(self, queries, subjects, evalue=0.001, max_rlen=0, command='blastn', **kwargs):
        queries_len = len(queries)
        subjects_len = len(subjects)
        results = [[None for _s in subjects] for _q in queries]
        pairs = list(itertools.product(xrange(queries_len), xrange(subjects_len)))
        @MultiprocessingBase.data_mapper
        @shelf_result
        def _worker(qs, queries, subjects):
            query = queries[qs[0]]
            subject = subjects[qs[1]]
            if query is None or subject is None: return None
            return LocalBlast.r2r_blast(query, subject, evalue, max_rlen, command, **kwargs)
        @MultiprocessingBase.results_assembler
        def _assembler(index, hsps, results, pairs):
            qs = pairs[index]
            results[qs[0]][qs[1]] = hsps
        work = self.Work()
        work.prepare_jobs(_worker, pairs, None, queries, subjects)
        work.set_assembler(_assembler, results, pairs)
        self.start_work(work)
        if not self.wait(work): return None
        return results
        
    @staticmethod
    @MultiprocessingBase.data_mapper
    @shelf_result
    def _translate_cds(fi, seq, table):
        f = seq.features[fi]
        sseq = f.extract(seq)
        try: 
            tseq = SeqRecord(sseq.seq.translate(table),
                             id=seq.id, name=seq.name,
                             description=seq.description)
            pf = SeqFeature(FeatureLocation(0, len(tseq)), 
                            id=f.id, type='CDS',
                            qualifiers=f.qualifiers)
            tseq.features.append(pf)
        except:
            print 'Unable to translate %s' % sseq 
            return None
        return tseq 
    
    @staticmethod
    @MultiprocessingBase.results_assembler
    def _translation_assembler(index, tname, translations):
        if tname and os.path.isfile(tname):
            with roDict(tname) as db:
                translations[index] = db['result']
            cleanup_file(tname)
    
    @staticmethod
    def percent_identity(hsp):
        return float(hsp.identities)/hsp.align_length
    
    @staticmethod
    @MultiprocessingBase.data_mapper
    def _find_features_by_hsps(qs, translations, blast_results):
        hsps_name = blast_results[qs[0]][qs[1]]
        if not hsps_name: return None
        with roDict(hsps_name) as db:
            hsps = db['result']
        if not hsps:
            print 'no hsps in %s: %s' % (hsps_name, hsps)  
            return None
        hsp = min(hsps, key=lambda h: h.expect)
        for f in translations[qs[1]].features:
            if hsp.sbjct_start in f or hsp.sbjct_end in f:
                return f, hsp.align_length, LocalBlast.percent_identity(hsp), hsp.expect
        return None
        
    @staticmethod
    def _features_of_type_i(record, ftype):
        return [i for i, f in enumerate(record.features) if f.type == ftype]
        
    @staticmethod
    def _get_features(ri, records):
        rec = records[ri]
        features = LocalBlast._features_of_type_i(rec, 'CDS')
        if not features:
            features = LocalBlast._features_of_type_i(rec, 'gene')
            if not features:
                print 'No gene/CDS features found in:'
                print rec.id, rec.description
                return None
        return features
    
    def _get_foi(self, records, foi):
        for q in foi: foi[q] = re.compile(foi[q])
        num_records = len(records)
        def _get_foi(ri, records, foi):
            r = records[ri]
            features = []
            for f in r.features:
                for q in foi:
                    qv = f.qualifiers.get(q)
                    if qv and foi[q].search(qv[0]):
                        features.append(f)
            return features
        features = self.parallelize_work(1, _get_foi, range(num_records), records, foi)
        return features
    
    _pgap = 'X'*20
    
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
        @param features_of_interest: dictionary of the form 
        {qualifier_name : qualifier_value}
        to mark features denoting known clusters that should be analyzed one 
        against the other
        @return: list of pairs (CDS_index, (blast_result1, blast_result2, ...)) 
        where CDS_index is an index of the gene/CDS feature in the reference.features list 
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
            features = self.parallelize_work(1, self._get_features, 
                                             range(num_records), 
                                             all_records)
            if not features or not features[0]:
                print 'Reference sequence does not contain annotated genes:'
                print reference.id, reference.description
                return None
            if len([f for f in features if f]) < 2:
                print 'Subject sequences do not contain annotated genes'
                return None
            #add gene ids
            for ri, genes in enumerate(features):
                if not genes: continue
                r = all_records[ri]
                for gene_id, gi in enumerate(genes):
                    r.features[gi].qualifiers['gene_id'] = gene_id
        #get features of interest if requested
        if features_of_interest:
            with user_message('Searching for features of interest...'):
                foi = self._get_foi(all_records, features_of_interest)
        else: foi = None
        #translate features to proteins
        with user_message('Translating genes found in the reference and subjects...'):
            translations = [None]*num_records
            for i, (f, rec) in enumerate(zip(features, all_records)):
                if not f: continue  
                translation = [None]*len(f)
                work = self.Work()
                work.prepare_jobs(self._translate_cds, f, None, rec, table)
                work.set_assembler(self._translation_assembler, translation)
                self.start_work(work)
                if not self.wait(work): return None
                if i > 0:
                    for t in translation:
                        if not t: continue
                        if translations[i]:
                            translations[i] += self._pgap+t
                        else: translations[i] = t
                else: translations[i] = translation
        #blast features against subjects
        with user_message('Performing local blast of every translated gene in the reference against every translated subject...'):
            stranslations = translations[1:]
            blast_results = self.r2r_blast_batch(translations[0], stranslations, evalue, 
                                                 max_rlen, command='blastp', task='blastp')
            if not blast_results: return None
        #parse blast results
        with user_message('Searching for genes in subjects that overlap with top blast hits...'):
            pairs = list(itertools.product(xrange(len(translations[0])), xrange(len(stranslations))))
            work = self.Work()
            work.prepare_jobs(self._find_features_by_hsps, pairs,
                              None, stranslations, blast_results)
            @MultiprocessingBase.results_assembler
            def _assembler(index, result, blast_results, pairs):
                qs = pairs[index]
                blast_results[qs[0]][qs[1]] = result
            work.set_assembler(_assembler, blast_results, pairs)
            self.start_work(work)
            if not self.wait(work): return None
        return zip((reference.features[f] for f in features[0]), blast_results)
    
    @staticmethod
    def g2g_to_csv(filename, g2g_results):
        locus_tag = lambda g, f, i: f.qualifiers.get('locus_tag', ['%s_%s' % (g.name, i+1)])[0]
        product = lambda f: f.qualifiers.get('product', [''])[0]
        with open(filename, 'wb') as out:
            writer = csv.writer(out)
            header = [g1.name+'_gene', 'locus_tag', g1.name+'_desc']
            for g in genomes:
                header += [g.name+'_gene', 'locus_tag', g.name+'_desc', g.name+'_coverage', g.name+'_percent', g.name+'_evalue']
            writer.writerow(header)
            for i, (feature, res) in enumerate(results):
                if not res: continue
                row = [i+1, locus_tag(g1, feature, i), product(feature)]
                for sres in res:
                    if sres is None:
                        row += ['', '', '', '', '', '']
                    else:
                        f, align_length, percent_identity, evalue = sres
                        gid = f.qualifiers.get('gene_id', 'N/A')
                        row += [gid, locus_tag(g, f, gid), product(f), 
                                '%d/%d' % (align_length, len(f)), percent_identity*100, evalue]
                writer.writerow(row)
#end class
        
#tests
import signal
from time import sleep
from DegenPrimer import tmpStorage

_pid = -1
abort_event = None
def sig_handler(signal, frame):
    if _pid != os.getpid(): return
    print('\nAborting. This may take some time '
          'as not all operations could be stopped immediately.\n')
    abort_event.set(); sleep(0.1)
    tmpStorage.clean_tmp_files()
#end def

if __name__ == '__main__':
    from multiprocessing import Event
    from BioUtils import CommonTools
    _pid = os.getpid()
    #setup signal handler
    signal.signal(signal.SIGINT,  sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)
    signal.signal(signal.SIGQUIT, sig_handler)
    
    abort_event = Event()
    lb = LocalBlast(abort_event)
    
    with CommonTools.user_message('Loading genomes...'):
        genomes_dir = u'/home/allis/Dropbox/Science/Микра/Thermococcus/sequence/GenBank/Thermococcus'
        genome_names = ['Thermococcus_barophilus_Ch5-complete.gb', 'Thermococcus_onnurineus_NA1-complete-genome.gb'] 
        genomes = CommonTools.load_files(abort_event, [os.path.join(genomes_dir, f) for f in genome_names], 'gb') 
    
    g1 = genomes[0]
    genomes = genomes[1:]
    
    @shelf_result
    def g2g2shelf():
        return lb.g2g_blastp(g1, genomes, 11, features_of_interest={'ugene_name': 'FC-full'})
        
    g2g_res = ''#'/tmp/DP-PCR-cz4Elt'
    if not os.path.isfile(g2g_res):
        g2g_res = g2g2shelf()
        print g2g_res
    if g2g_res:
        with roDict(g2g_res) as db:
            results = db['result'] 
        if results:
            lb.g2g_to_csv('g2g_test.csv', results)
    
    print 'Done.'
