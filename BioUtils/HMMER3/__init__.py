# coding=utf-8

'''
Created on Dec 25, 2015

@author: Allis Tauri <allista@gmail.com>
'''

import sys, os
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
from Bio.SearchIO.HmmerIO import Hmmer3TextParser
from Bio.SeqFeature import SeqFeature, FeatureLocation

from BioUtils.Tools.Multiprocessing import MultiprocessingBase, cpu_count
from BioUtils.Tools.tmpStorage import shelf_result, roDict

from BioUtils.Tools import retry, mktmp_name
from BioUtils.Tools.Output import user_message, Progress, ProgressCounter
from BioUtils.SeqUtils import mktmp_fasta, cat_records, Translator, get_indexes_of_genes
from BioUtils.HMMER3.Applications import HMMSearchCommandline

class BatchHmmer(MultiprocessingBase):
    
    def __init__(self, abort_event):
        super(BatchHmmer, self).__init__(abort_event)
    
    @staticmethod
    def hmmsearch_recs(hmm, recs, **kwargs):
        recfile = mktmp_fasta(recs)
        hmm_out = mktmp_name('.hmm.txt')
        try:
            cline = HMMSearchCommandline(hmmfile=hmm, seqdb=recfile,
                                         o=hmm_out, cpu=cpu_count, seed=0,
                                         **kwargs)
            stdout, stderr = cline()
            print stdout
            if stderr: 
                sys.stderr.write(stderr)
                sys.stderr.flush()
            #parse hmmsearch results
            with open(hmm_out) as inp:
                return list(Hmmer3TextParser(inp))
        except Exception, e:
            print 'Error while running hmmsearch.'
            print e
            return None
        finally:
            os.unlink(recfile)
            os.unlink(hmm_out)
    
    def hmmsearch_genome(self, hmm, genome, table='Standard', decorate=False, **kwargs):
        #get genes
        genes = get_indexes_of_genes(genome)
        if not genes: return None
        for gene_id, gi in enumerate(genes):
            genome.features[gi].qualifiers['feature_id'] = gi
            genome.features[gi].qualifiers['gene_id'] = gene_id
        #translate genes
        with user_message('Translating genes/CDS of %s' % genome.description, '\n'):
            translator = Translator(self._abort_event)
            translation = translator.translate(genome, genes, table)
        if not translation: return None
        with user_message('Performing hmm search.'):
            results = self.hmmsearch_recs(hmm, translation)
        if not results: return None
        with user_message('Parsing search results...'):
            #get hit_ids of hmm matches
            hits = dict()
            for result in results:
                for hit in result.iterhits():
                    hits[hit.id] = hit
            #get indexes of features where hmm hit
            hit_features = dict()
            for t in translation:
                if t.id in hits:
                    fid = t.features[0].qualifiers.get('feature_id')
                    if fid is None: continue
                    hit_features[fid] = hits[t.id], t
        #decorate genome
        if decorate:
            with user_message('Adding results as annotations...'):
                hmm_name = os.path.basename(hmm)
                for f in hit_features:
                    feature = genome.features[f] 
                    for hsp in hit_features[f][0]:
                        if feature.strand == 1:
                            hmm_location = FeatureLocation(feature.location.start+hsp.hit_start*3, 
                                                           feature.location.start+hsp.hit_end*3, 
                                                           feature.strand)
                        else:
                            hmm_location = FeatureLocation(feature.location.end-hsp.hit_end*3, 
                                                           feature.location.end-hsp.hit_start*3, 
                                                           feature.strand)
                        hmm_feature = SeqFeature(hmm_location, type='misc_feature')
                        hmm_feature.qualifiers['hmm_model'] = hmm_name
                        hmm_feature.qualifiers['bitscore'] = hsp.bitscore
                        hmm_feature.qualifiers['psi_evalue'] = hsp.psi_evalue
                        hmm_feature.qualifiers['evalue_cond'] = hsp.evalue_cond
                        hmm_feature.qualifiers['acc_average'] = hsp.acc_avg
                        hmm_feature.qualifiers['bias'] = hsp.bias
                        genome.features.append(hmm_feature)
        print 'Done.\n'
        return hit_features 
    
#tests
import signal
from time import sleep

from BioUtils.NCBI import BlastCLI
from reportlab.lib import colors

from BioUtils.Tools.tmpStorage import roDict, clean_tmp_files, shelf_result

_pid = -1
abort_event = None
def sig_handler(signal, frame):
    if _pid != os.getpid(): return
    print('\nAborting. This may take some time '
          'as not all operations could be stopped immediately.\n')
    abort_event.set(); sleep(0.1)
    clean_tmp_files()
#end def

if __name__ == '__main__':
    from multiprocessing import Event
    from BioUtils.Tools.Output import user_message
    from BioUtils.SeqUtils import load_files, load_dir
    _pid = os.getpid()
    #setup signal handler
    signal.signal(signal.SIGINT,  sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)
    signal.signal(signal.SIGQUIT, sig_handler)
    
    abort_event = Event()
    with user_message('Loading genomes...', '\n'):
        genomes_dir = u'/home/allis/Documents/INMI/Aerobic-CODH/genomes/'
        genome_names = ['Thermococcus_barophilus_Ch5-complete.gb', 
                        'Thermococcus_onnurineus_NA1-complete-genome.gb',
                        'Thermococcus_sp._ES1.gb',
                        'Thermococcus-DS1-preliminary.gb'] 
        genomes = load_dir(abort_event, genomes_dir, 'gb', r'.*\.gb')
        if not genomes: sys.exit(1)
#        load_files(abort_event, [os.path.join(genomes_dir, f) for f in genome_names], 'gb') 
    
    hmm = u'/home/allis/Documents/INMI/Aerobic-CODH/COX-EC/COX-EC_1.2.99.2_CoxL.hmm'
    
    hmmer = BatchHmmer(abort_event)
    
    for g in genomes:
        results = hmmer.hmmsearch_genome(hmm, g, table=11, decorate=True)
        if results: 
            SeqIO.write(g, '%s.gb' % g.name, 'gb')
            print '='*80
            print g.name, g.description
            for fi in results:
                print results[fi][0]
                print results[fi][1]
                print
    print '\nDone'