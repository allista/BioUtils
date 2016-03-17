# coding=utf-8

'''
Created on Mar 17, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os
import re
import csv
import itertools

from random import shuffle

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast import Applications as BPApps

from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from BioUtils.Tools.tmpStorage import shelf_result, roDict

from BioUtils.Tools.Misc import mktmp_name, safe_unlink
from BioUtils.Tools.Output import user_message, Progress, ProgressCounter
from BioUtils.SeqUtils import mktmp_fasta, cat_records, Translator

from .Applications import BlastDBcmdCommandline
from .BlastBase import _BlastBase

class BlastCLI(MultiprocessingBase, _BlastBase):

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
        _BlastBase.__init__(self)
    
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
    
    def fetch_results(self, results, db, what='record'):
        '''Fetch records that were found by BLAST search from the database
        @param results: an iterable of Bio.Blast.Record.Blast objects (aka BlastRecords)
        @param db: name of the BLAST database to get records from
        @param what: string, one of "record", "alignment" or "hsp"
        @return: list of SeqRecord objects'''
        #parse results into queries
        queries = []
        for record in results:
            for alignment in record.alignments:
                q = self.Query(alignment, what)
                if not q: continue
                queries.append(q)
        if not queries: return None
        #fetch records from database
        out_file = mktmp_name('.fasta')
        entry_batch = mktmp_name('.qry')
        qstr = '\n'.join(str(q) for q in queries)
        with open(entry_batch, 'w') as out: out.write(qstr)
        try:
            cline = BlastDBcmdCommandline(db=db, out=out_file, 
                                          entry_batch=entry_batch)
            out, err = cline()
            if err:
                print '\nError in blastdbcmd call'
                print cline
                print out
                print err
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
                print ('\nReference sequence does not contain annotated genes:\n%s %s' 
                       % (reference.id, reference.description))
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
