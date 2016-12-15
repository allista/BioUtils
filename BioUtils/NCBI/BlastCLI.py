# coding=utf-8

'''
Created on Mar 17, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os
import re
import csv
import itertools
import tempfile
import shutil

from itertools import chain
from random import shuffle

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast import Applications as BPApps

from BioUtils.Tools.Debug import estr
from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from BioUtils.Tools.tmpStorage import shelf_result, roDict, from_shelf

from BioUtils.Tools.Text import random_text
from BioUtils.Tools.Misc import mktmp_name, safe_unlink, run_cline
from BioUtils.Tools.Output import user_message, Progress, ProgressCounter
from BioUtils.SeqUtils import mktmp_fasta, cat_records, Translator, pretty_rec_name
from BioUtils.Annotation import HSPAnnotator

from .Applications import BlastDBcmdCommandline, FormatDBCommandline
from .BlastBase import _BlastBase

class BlastCLI(MultiprocessingBase, _BlastBase, HSPAnnotator):

    _clines = {
            'blastp'     : BPApps.NcbiblastpCommandline, 
            # compares an amino acid query sequence against a protein sequence database
            
            'blastn'     : BPApps.NcbiblastnCommandline, 
            # compares a nucleotide query sequence against a nucleotide sequence database
            
            'blastx'     : BPApps.NcbiblastxCommandline, 
            # compares a nucleotide query sequence translated in all reading frames against a protein sequence database
            
            'tblastn'    : BPApps.NcbitblastnCommandline, 
            # compares a protein query sequence against a nucleotide sequence database dynamically translated in all reading frames
            
            'tblastx'    : BPApps.NcbitblastxCommandline, 
            # compares the six-frame translations of a nucleotide query sequence against the six-frame translations of a nucleotide sequence database. Please note that tblastx program cannot be used with the nr database on the BLAST Web page.
            
            'psiblast'   : BPApps.NcbipsiblastCommandline, # Position-Specific Initiated BLAST
            'rpsblast'   : BPApps.NcbirpsblastCommandline, # Reverse Position Specific BLAST
            'rpstblastn' : BPApps.NcbirpstblastnCommandline, # Translated Reverse Position Specific BLAST
            'deltablast' : BPApps.NcbideltablastCommandline # Protein-Protein domain enhanced lookup time accelerated blast
            }

    def __init__(self, abort_event):
        super(BlastCLI, self).__init__(abort_event)
        _BlastBase.__init__(self, abort_event)
        HSPAnnotator.__init__(self)
    
    @classmethod
    def blast(cls, command, **kwargs):
        '''Generic wrapper for blast commandline programs'''
        results_file = kwargs.pop('save_results_to', None)
        parse_results = kwargs.pop('parse_results', False)
        if results_file: bout = results_file
        else:
            parse_results = True
            bout = mktmp_name('.xml')
        try:
            cline = cls._clines[command]
            cmd = cline(outfmt=5, out=bout, **kwargs)
            out, err = cmd()
            if err:
                print ('\nError message from %s:\n%s\n%s\n%s'
                       % (command, cmd, out, err))
            if parse_results:
                with open(bout) as inp:
                    results = list(NCBIXML.parse(inp))
                return results
        except Exception, e:
            print '\nError while performing local %s:\n%s' % (command, estr(e))
            return None
        finally:
            if not results_file: os.unlink(bout)
        return None
    
    @classmethod
    def blast_seq(cls, query, db, evalue=0.001, command='blastn', **kwargs):
        '''Perform local blast of a SeqRecord against a database'''
        qfile = mktmp_fasta(query)
        results = cls.blast(command, query=qfile, db=db, 
                            evalue=evalue, **kwargs)
        os.unlink(qfile)
        return results
    
    @classmethod
    def format_tmp_db(cls, sequences, nucleotide=True):
        '''Create a temporary local blast database
        @param sequences: SeqRecord object to populate the database with
        @param nucleotide: if the sequences are nucleotide (True) or protein (False)
        @return: database name that includes includes its path'''
        basename = random_text(8)
        dbdir = tempfile.mkdtemp('_blastDB')
        sfile = mktmp_fasta(sequences)
        if run_cline(FormatDBCommandline(input=sfile, 
                                         protein='F' if nucleotide else 'T', 
                                         name=basename),
                     cwd=dbdir): 
            return os.path.join(dbdir, basename)
        else: shutil.rmtree(dbdir, ignore_errors=True)
        return None

    @classmethod
    def delete_tmp_db(cls, dbpath):
        dbdir = os.path.dirname(dbpath)
        shutil.rmtree(dbdir, ignore_errors=True)
    
    @staticmethod
    def base_sid(s): return s.id.split(':')[0]
    
    @classmethod
    def unique_seqs(cls, seqs):
        unique = {}
        for s in seqs:
            sid = cls.base_sid(s)
            if sid in unique:
                #FIXME: need to merge sequences properly, instead of replacing 
                if len(unique[sid]) < len(s): unique[sid] = s 
            else: unique[sid] = s
        return unique.values()
    
    def ring_blast(self, query, db='nr', evalue=0.001, blast_filter=None, depth=1, command='blastn', **kwargs):
        '''Perform a blast search with the given query to obtain the core set of hits.
        Make another search with each hit as a query.
        If results of the second search contain new hits,
        check if these are reciprocal by yet another search with them
        and checking that results contain hits from the core set and if they are,
        add the to the final set.
        '''
        if isinstance(query, SeqRecord): query = [query]
        def blast_filter_fetch(seqs):
            @MultiprocessingBase.data_mapper
            @shelf_result
            def worker(s):
                r = self.blast_seq(s, db, evalue, command, **kwargs)
                if r and blast_filter: blast_filter(r)
                if r: return self.fetch_results(r, db, what='alignment')
                return None
            results = []
            total = len(seqs)
            prg = ProgressCounter('Performing blast search for %d sequences:' % total, total)
            @MultiprocessingBase.results_assembler
            def assembler(i, res):
                if res: results.append(res)
                prg.count()
            with prg:
                if not self.parallelize2(1, worker, assembler, seqs): return None
                return results
        
        with user_message('RingBlast: building a core set of sequences.', '\n'):
            core_seqs = blast_filter_fetch(query)
            if not core_seqs: return None
            core_seqs = self.unique_seqs(chain.from_iterable(from_shelf(r) for r in core_seqs))
            extended_set = dict((self.base_sid(s), s) for s in core_seqs)
            if depth <= 0: return core_seqs
            core_db = self.format_tmp_db(core_seqs, command.endswith('n'))
            
        def check_sequences(seqs, next_to_process):
            total = len(seqs)
            prg = ProgressCounter('RingBlast: checking %d new sequences:' % total, total)
            @MultiprocessingBase.data_mapper
            def worker(seq):
                res = self.blast_seq(seq, core_db, 100, command)
                if res and blast_filter: blast_filter(res)
                return bool(res), seq
            @MultiprocessingBase.results_assembler
            def assembler(i, res):
                prg.count()
                if not res[0]: return 
                seq = res[1]
                extended_set[self.base_sid(seq)] = seq
                next_to_process.append(seq)
            with prg: return self.parallelize2(1, worker, assembler, seqs)
            
        def process_sequences(seqs, _depth):
            if _depth == 0: return
            with user_message('RingBlast: processing %d sequences of the %d ring.' 
                              % (len(seqs), depth-_depth+1), '\n'): 
                next_ring = blast_filter_fetch(seqs)
                if not next_ring: return
                to_check = []
                next_to_process = []
                for n in next_ring:
                    next_seqs = from_shelf(n)
                    if not next_seqs: continue 
                    for ns in next_seqs:
                        sid = self.base_sid(ns)
                        if sid in extended_set:
                            #FIXME: need to merge sequences properly, instead of replacing 
                            if len(extended_set[sid]) < len(ns):
                                extended_set[sid] = ns 
                        else: to_check.append(ns)
            if not to_check or not check_sequences(to_check, next_to_process): return
            if next_to_process: process_sequences(next_to_process, _depth-1)
            
        process_sequences(core_seqs, depth)
        return extended_set.values()

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
    def fetch_queries(cls, queries, db):
        out_file = mktmp_name('.fasta')
        entry_batch = mktmp_name('.qry')
        qstr = '\n'.join(str(q) for q in queries)
        with open(entry_batch, 'w') as out: out.write(qstr)
        try:
            cline = BlastDBcmdCommandline(db=db, out=out_file,
                                          entry_batch=entry_batch)
            out, err = cline()
            if err:
                print '\nError in blastdbcmd call:\n%s\n%s\n%s' % (cline, out, err)
                return None
            #parse results
            records = list(SeqIO.parse(out_file, 'fasta'))
            return records
        except Exception, e:
            print e
            return None
        finally:
            safe_unlink(entry_batch)
            safe_unlink(out_file)

    @classmethod
    def fetch_results(cls, results, db, what='record'):
        '''Fetch records that were found by BLAST search from the database
        @param results: an iterable of Bio.Blast.Record.Blast objects (aka BlastRecords)
        @param db: name of the BLAST database to get records from
        @param what: string, one of "record", "alignment" or "hsp"
        @return: list of SeqRecord objects'''
        #parse results into queries
        queries = []
        for record in results:
            for alignment in record.alignments:
                q = cls.Query(alignment, what)
                if not q: continue
                queries.append(q)
        if not queries: return None
        #fetch records from database
        return cls.fetch_queries(queries, db)

    def _s2s_blast_batch(self, queries, subjects, subject_locs=None, evalue=0.001, command='blastn', **kwargs):
        queries_len = len(queries)
        subjects_len = len(subjects)
        results = [[None for _s in subjects] for _q in queries]
        pairs = list(itertools.product(xrange(queries_len), xrange(subjects_len)))
        ignore_none_locs = kwargs.pop('ignore_none_locs', False)
        shuffle(pairs)
        @MultiprocessingBase.data_mapper
        @shelf_result
        def worker(qs, queries, subjects, subject_locs):
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
        def assembler(index, hsps, results, pairs, prg):
            qs = pairs[index]
            results[qs[0]][qs[1]] = hsps
            prg.count()
        with ProgressCounter('Performing multiple %s searches:'%command, len(pairs)) as prg:
            work = self.Work()
            work.start_work(worker, pairs, None, queries, subjects, subject_locs)
            work.assemble(assembler, results, pairs, prg)
            if not work.wait(): return None
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
                print 'No gene/CDS features found in:\n%s %s' % (rec.id, rec.description)
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
                print ('\nReference sequence does not contain annotated _genes:\n%s %s'
                       % (reference.id, reference.description))
                return None
            if len([f for f in features if f]) < 2:
                print '\nSubject sequences do not contain annotated _genes'
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
        with Progress('Translating _genes found in the reference and subjects...', num_records) as prg:
            translator = Translator(self._abort_event)
            translations = [None]*num_records
            foi_translations = [[None]*num_records for _f in fois]
            for i, (f, rec) in enumerate(zip(features, all_records)):
                if not f:
                    prg.step(i) 
                    continue
                translation = translator.translate_features(rec, f, table)
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
        with ProgressCounter('Searching for _genes in subjects that overlap with top blast hits...', len(pairs)) as prg:
            work = self.Work()
            work.start_work(self._find_features_by_hsps, pairs,
                            None, stranslations, blast_results)
            @MultiprocessingBase.results_assembler
            def assembler(index, result, blast_results, pairs, prg):
                qs = pairs[index]
                blast_results[qs[0]][qs[1]] = result
                prg.count()
            work.assemble(assembler, blast_results, pairs, prg)
            if not work.wait(): return None
        return zip((reference.features[f] for f in features[0]), blast_results)
    
    @staticmethod
    def g2g_to_csv(filename, reference, subjects, g2g_results):
        locus_tag = lambda g, f, i: f.qualifiers.get('locus_tag', ['%s_%s' % (g.name, i+1)])[0]
        product = lambda f: f.qualifiers.get('product', [''])[0]
        with open(filename, 'wb') as out:
            writer = csv.writer(out)
            header = [reference.name+'gene', 'locus_tag', reference.name+'_desc']
            for g in subjects:
                header += [g.name+'gene', 'locus_tag', g.name+'_desc', g.name+'_coverage', g.name+'_percent', g.name+'_evalue']
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

    @staticmethod
    def _get_homologues(query, min_identity, db, command):
        homologues = set()
        qlen = float(len(query))
        results = BlastCLI.blast_seq(query, db, evalue=10, command=command)
        if not results: return homologues
        for r in results:
            for ali in r.alignments:
                identities = sum(hsp.identities for hsp in ali.hsps)
                if identities / qlen > min_identity:
                    rid = ali.hit_def.split()[0]
                    if rid: homologues.add(rid)
        try: homologues.remove(query.id)
        except KeyError: pass
        return homologues

    @staticmethod
    def _filter_homologues(get_all_homologues, seqs, min_identity, keep_ids=None, nucleotide=False):
        print 'Filtering out close homologues. This will take a wile:'
        command = 'blastn' if nucleotide else 'blastp'
        dbname = ''
        try:
            with user_message('Formatting blast DB', '\n'):
                dbname = BlastCLI.format_tmp_db(seqs, nucleotide)
                if not dbname:
                    print 'Unable to make temporary BLAST database.'
                    return None
            with ProgressCounter('Searching for homologues using local blastp...', len(seqs)) as prg:
                homologues = get_all_homologues(seqs, min_identity, dbname, command, prg)
        except Exception as e:
            print '%s\n' % str(e)
            return None
        finally:
            if dbname:
                shutil.rmtree(os.path.dirname(dbname), ignore_errors=True)
        if not homologues: return seqs
        with user_message('Removing all homologs from each group except the first one...'):
            remove = set()
            if keep_ids: keep_ids = set(keep_ids)
            for seq in seqs:
                if seq.id in remove: continue
                h = homologues.pop(seq.id, set())
                if h:
                    if keep_ids:
                        nhoms = len(h)
                        h -= keep_ids
                        if nhoms != len(h) and seq.id not in keep_ids:
                            h.add(seq.id)
                    remove.update(h)
            return [seq for seq in seqs if seq.id not in remove]

    def filter_homologues(self, seqs, min_identity, keep_ids=None, nucleotide=False):
        if not seqs: return False
        @MultiprocessingBase.data_mapper
        def _worker(qi, queries, db, command):
            query = queries[qi]
            homologues = BlastCLI._get_homologues(query, min_identity, db, command=command)
            return query.id, homologues
        @MultiprocessingBase.results_assembler
        def _assembler(qi, result, homologues, prg):
            homologues[result[0]] = result[1]
            prg.count()
        def get_all_homologues(seqs, min_identity, dbname, command, prg):
            homologs = {}
            work = self.Work()
            work.start_work(_worker, range(len(seqs)), None, seqs, dbname, command)
            work.assemble(_assembler, homologs, prg)
            if not work.wait(): return {}
            return homologs
        return self._filter_homologues(get_all_homologues, seqs, min_identity, keep_ids, nucleotide)

    def filter_homologues_single(self, seqs, min_identity, keep_ids=None, nucleotide=False):
        if not seqs: return False
        def get_all_homologues(seqs, min_identity, dbname, command, prg):
            homologues = {}
            for seq in seqs:
                if self.aborted(): return {}
                seq_homologues = set()
                qlen = float(len(seq))
                results = BlastCLI.blast_seq(seq, dbname, evalue=10, command=command)
                if not results:
                    prg.count()
                    continue
                for r in results:
                    for ali in r.alignments:
                        identities = sum(hsp.identities for hsp in ali.hsps)
                        if identities / qlen > min_identity:
                            rid = ali.hit_def.split()[0]
                            if rid: seq_homologues.add(rid)
                try: seq_homologues.remove(seq.id)
                except KeyError: pass
                homologues[seq.id] = seq_homologues
                prg.count()
            return homologues
        return BlastCLI._filter_homologues(get_all_homologues, seqs, min_identity, keep_ids, nucleotide)

    annotation_type = 'blast'

    def hsp_score(self, hsp):
        return hsp.identities/float(hsp.align_length)*100

    def hsp2feature(self, name, group, location, hsp):
        feature = super(BlastCLI, self).hsp2feature(name, group, location, hsp)
        feature.qualifiers['bitscore'] = hsp.bits
        feature.qualifiers['evalue'] = hsp.expect
        feature.qualifiers['gaps'] = hsp.gaps
        feature.qualifiers['identities'] = hsp.identities
        feature.qualifiers['align_length'] = hsp.align_length
        feature.qualifiers['percent'] = self.hsp_score(hsp)
        feature.qualifiers['query_start'] = hsp.query_start
        feature.qualifiers['query_end'] = hsp.query_end
        return feature

    @staticmethod
    def add_program(feature, program):
        feature.qualifiers['program'] = program

    def blastn_annotate(self, tag_sequences, subject_record, min_identity, evalue=0.001, **kwargs):
        results = self.s2s_blast_batch(tag_sequences, [subject_record], evalue=evalue, command='blastn', **kwargs)
        if results is None: return False
        with user_message('Adding results as annotations...'):
            annotated = False
            for i, tag in enumerate(tag_sequences):
                if not results[i]: continue
                record = results[i][0]
                if not record: continue
                tag_name = pretty_rec_name(tag)
                if tag_name != tag.id:
                    tag_name += ' (%s)' % tag.id
                for hit in record:
                    for ali in hit.alignments:
                        for hsp in ali.hsps:
                            if hsp.identities / float(hsp.align_length) < min_identity: continue
                            strand = 1 if hsp.sbjct_start < hsp.sbjct_end else -1
                            if strand == 1:
                                location = FeatureLocation(hsp.sbjct_start-1,
                                                           hsp.sbjct_end,
                                                           strand)
                            else:
                                location = FeatureLocation(hsp.sbjct_end-1,
                                                           hsp.sbjct_start,
                                                           strand)
                            feature = self.hsp2feature(tag_name,'blastn_annotations', location, hsp)
                            self.add_program(feature, 'blastn')
                            subject_record.features.append(feature)
                            annotated = True
        return annotated

    def blastp_annotate(self, tag_sequences, subject_record, min_identity, evalue=0.001, table=11, **kwargs):
        # translate subject in six frames
        with user_message('Translating whole genome in 6 reading frames', '\n'):
            translator = Translator(self._abort_event)
            translation = translator.translate_six_frames(subject_record, table)
        if not translation: return False
        results = self.s2s_blast_batch(tag_sequences, translation, evalue=evalue, command='blastp', **kwargs)
        if results is None: return False
        with user_message('Adding results as annotations...'):
            annotated = False
            subj_len = len(subject_record)
            for i, tag in enumerate(tag_sequences):
                if not results[i]: continue
                tag_name = pretty_rec_name(tag)
                if tag_name != tag.id:
                    tag_name += ' (%s)' % tag.id
                for frame, record in enumerate(results[i]):
                    if not record: continue
                    frec = translation[frame]
                    start = frec.annotations['start']
                    strand = frec.annotations['strand']
                    for hit in record:
                        for ali in hit.alignments:
                            for hsp in ali.hsps:
                                if hsp.identities / float(hsp.align_length) < min_identity: continue
                                if strand == 1:
                                    location = FeatureLocation(start+(hsp.sbjct_start-1)*3,
                                                               start+hsp.sbjct_end*3,
                                                               strand)
                                else:
                                    location = FeatureLocation(subj_len-start-hsp.sbjct_end*3,
                                                               subj_len-start-hsp.sbjct_start*3,
                                                               strand)
                                feature = self.hsp2feature(tag_name, 'blastp_annotations', location, hsp)
                                self.add_program(feature, 'blastp')
                                subject_record.features.append(feature)
                                annotated = True
        return annotated

#end class
