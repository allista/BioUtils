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
    
    def r2r_blast_batch(self, queries, subjects, subject_locs=None, evalue=0.001, max_rlen=0, command='blastn', **kwargs):
        results = self._r2r_blast_batch(queries, subjects, subject_locs, evalue, max_rlen, command, **kwargs)
        if not results: return None
        for qi in xrange(len(results)):
            for si in xrange(len(results[qi])):
                hsps_name = results[qi][si]
                if not hsps_name: continue
                with roDict(hsps_name) as db:
                    results[qi][si] = db['result']
        return results
    
    def _r2r_blast_batch(self, queries, subjects, subject_locs=None, evalue=0.001, max_rlen=0, command='blastn', **kwargs):
        queries_len = len(queries)
        subjects_len = len(subjects)
        results = [[None for _s in subjects] for _q in queries]
        pairs = list(itertools.product(xrange(queries_len), xrange(subjects_len)))
        @MultiprocessingBase.data_mapper
        @shelf_result
        def _worker(qs, queries, subjects, subject_locs):
            query = queries[qs[0]]
            subject = subjects[qs[1]]
            if query is None or subject is None: return None
            if subject_locs and subject_locs[qs[1]]: 
                loc = tuple(subject_locs[qs[1]])
                kwargs['subject_loc'] = '%d-%d' % loc
            elif 'subject_loc' in kwargs: del kwargs['subject_loc']
            return LocalBlast.r2r_blast(query, subject, evalue, max_rlen, command, **kwargs)
        @MultiprocessingBase.results_assembler
        def _assembler(index, hsps, results, pairs):
            qs = pairs[index]
            results[qs[0]][qs[1]] = hsps
        work = self.Work()
        work.prepare_jobs(_worker, pairs, None, queries, subjects, subject_locs)
        work.set_assembler(_assembler, results, pairs)
        self.start_work(work)
        if not self.wait(work): return None
        return results
        
    @staticmethod
    @MultiprocessingBase.data_mapper
    @shelf_result
    def _translate_genes(fi, rec, table):
        f = rec.features[fi]
        srec = f.extract(rec)
        try: 
            tsec = srec.seq.translate(table)
            if tsec[-1] == '*': tsec = tsec[:-1]
            trec = SeqRecord(tsec,
                             id=rec.id, name=rec.name,
                             description=rec.description)
            pf = SeqFeature(FeatureLocation(0, len(trec)), 
                            id=f.id, type='CDS',
                            qualifiers=f.qualifiers)
            trec.features.append(pf)
        except:
            raise RuntimeError('Unable to translate: %s' % str(srec.seq))
        return trec 
    
    @staticmethod
    @MultiprocessingBase.results_assembler
    def _translation_assembler(index, tname, translations):
        if tname and os.path.isfile(tname):
            with roDict(tname) as db:
                translations[index] = db['result']
            cleanup_file(tname)
    
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
#        print '\n'+'='*80
#        for h in hsps:
#            if h is hsp:
#                print '*'*80
#                print h
#                print '*'*80
#            else: print h
#            print
#        print '-'*80
        for f in translations[qs[1]].features:
#            print f.qualifiers.get('locus_tag'), f.location, '%d-%d' % (hsp.sbjct_start, hsp.sbjct_end), \
#                  'MATCH' if hsp.sbjct_start-1 in f and hsp.sbjct_end-1 in f else ''
            if hsp.sbjct_start-1 in f and hsp.sbjct_end-1 in f:
                aln_len = abs(hsp.sbjct_start-hsp.sbjct_end)+1
                return f, aln_len, float(hsp.identities)/aln_len, hsp.expect
        print '\n'+'='*80
        print
        return None
        
    @staticmethod
    def _features_of_type_i(record, ftype):
        return [i for i, f in enumerate(record.features) if f.type == ftype]
        
    @staticmethod
    def _get_genes(rec):
        features = LocalBlast._features_of_type_i(rec, 'CDS')
        if not features:
            features = LocalBlast._features_of_type_i(rec, 'gene')
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
    
    @staticmethod
    def _cat_records(records, gap='X', glen=20):
        cat = None
        gap = gap*glen
        for t in records:
            if t: 
                if cat: cat += gap+t
                else: cat = t
        return cat
    
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
            features = self.parallelize_work(1, lambda ri, records: self._get_genes(records[ri]), 
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
                    r.features[gi].qualifiers['feature_id'] = gi
                    r.features[gi].qualifiers['gene_id'] = gene_id
        #get features of interest if requested
        fois = None
        if features_of_interest:
            with user_message('Searching for features of interest...'):
                fois = []
                for foi in features_of_interest:
                    fois.append(self._get_fois(all_records, foi))
        #translate features to proteins
        with user_message('Translating genes found in the reference and subjects...'):
            translations = [None]*num_records
            foi_translations = [[None]*num_records for _f in fois]
            for i, (f, rec) in enumerate(zip(features, all_records)):
                if not f: continue
                translation = [None]*len(f)
                work = self.Work()
                work.prepare_jobs(self._translate_genes, f, None, rec, table)
                work.set_assembler(self._translation_assembler, translation)
                self.start_work(work)
                if not self.wait(work): return None
                if i > 0: 
                    translations[i] = self._cat_records(translation)
                    if fois:
                        for ifoi, foi in enumerate(fois):
                            foi_loc = [0, 0]
                            for foi_var in foi[i]: 
                                for gid in foi_var:
                                    l = translations[i].features[gid].location
                                    foi_loc[0] = min(int(l.start)+1, foi_loc[0]) if foi_loc[0] > 0 else int(l.start)+1
                                    foi_loc[1] = max(int(l.end), foi_loc[1])
                            if foi_loc[0] >= 0: foi_translations[ifoi][i] = foi_loc 
                else: 
                    translations[i] = translation
                    if fois: 
                        for ifoi, foi in enumerate(fois):
                            foi_translations[ifoi][i] = [[translation[gid] for gid in foi_var] for foi_var in foi[i]]
#DEBUG
#        translations[1].name = 'TON'
#        ch5 = self._cat_records(translations[0])
#        ch5.name = 'TBCH5'
#        SeqIO.write([ch5, translations[1]], 'trans-test.gb', 'gb')
        #blast features against subjects
        with user_message('Performing local blast of every translated gene in the reference against every translated subject...'):
            stranslations = translations[1:]
            blast_results = self._r2r_blast_batch(translations[0], stranslations, None, evalue, 
                                                  max_rlen, command='blastp', task='blastp')
            if not blast_results: return None
            if fois: #redo blast for fois and replace the results
                for ifoi, foi in enumerate(foi_translations):
                    sfoi_locs = foi[1:]
                    for i, foi_var in enumerate(foi[0]):
                        foi_blast = self._r2r_blast_batch(foi_var, stranslations, sfoi_locs, evalue, 
                                                          max_rlen, command='blastp', task='blastp')
                        if not foi_blast: continue
                        for gi, gid in enumerate(fois[ifoi][0][i]):
                            blast_results[gid] = foi_blast[gi]
        #process blast results
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
import signal, sys
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
    
#    from DegenPrimer import MultiprocessingBase
#    MultiprocessingBase.cpu_count = 1
    abort_event = Event()
    lb = LocalBlast(abort_event)
    
    with CommonTools.user_message('Loading genomes...'):
        genomes_dir = u'/home/allis/Dropbox/Science/Микра/Thermococcus/sequence/GenBank/Thermococcus'
        genome_names = ['Thermococcus_barophilus_Ch5-complete.gb', 
                        'Thermococcus_onnurineus_NA1-complete-genome.gb',
                        'Thermococcus_sp._ES1.gb',
                        'Thermococcus-DS1-preliminary.gb'] 
        genomes = CommonTools.load_files(abort_event, [os.path.join(genomes_dir, f) for f in genome_names], 'gb') 
    
#    g1 = genomes[0][1460915:1521600]
#    genomes = [genomes[1][1408604:1478016]]
    
    g1 = genomes[0]
    genomes = genomes[1:]
    
    @shelf_result
    def g2g2shelf():
        return lb.g2g_blastp(g1, genomes, 11, features_of_interest=[{'ugene_name': 'FC-full'}, {'ugene_name': 'COC-full'}])
        
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
