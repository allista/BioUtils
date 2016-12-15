# coding=utf-8

'''
Created on Mar 16, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os

from BioUtils.Tools.Text import FilenameParser
from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.AlignmentUtils import AlignmentUtils
from BioUtils.SeqUtils import copy_attrs, SeqView
from BioUtils.Tools.Output import user_message
from BioUtils.PrimerFinder import PrimerFinder


class Main(MPMain):
    def _main(self):
        min_prod = 400
        silva_db = '/home/allis/Documents/INMI/SILVA-DB/SILVA_123_SSURef_Nr99_tax_silva.fasta'
        alifile = '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.aln.fasta'
        add_filename = FilenameParser.strip_ext(alifile)+'.with_additions.fasta'
        outgroups = ['Thermococcus_chitonophagus', 'SMTZ1-55', 'contig72135_1581_sunspring_meta']
        add = ['KF836721.1.1270','EU635905.1.1323']
        exclude = []#['Thermococcus_chitonophagus', 'SMTZ1-55', 'BA1-16S', 'contig72135_1581_sunspring_meta']
        #load alignment
        if os.path.isfile(add_filename): 
            alifile = add_filename
            add_filename = ''
        with user_message('Loadding initial alignment...', '\n'):
            orig_ali = AlignmentUtils.load_first(alifile)
            if not orig_ali: return 1
        #load homologs
        if add_filename:
            with user_message('Loadding additional sequences...', '\n'):
                add_seqs = []
                db = SeqView()
                if db.load(silva_db):
                    for sid in add:
                        seq = db.get(sid)
                        if seq: add_seqs.append(seq)
                        else: print '%s not found in %s' % (sid, silva_db)
            #realign data if needed
            if add_seqs:
                with user_message('Realigning data...', '\n'):
                    add_filename = FilenameParser.strip_ext(alifile)+'.with_additions.fasta'
                    AlignmentUtils.align(list(orig_ali)+add_seqs, add_filename)
                    orig_ali = AlignmentUtils.load_first(add_filename)
                    if not orig_ali: return 2
        #process the alignment
        ali = orig_ali.remove(*exclude).trim()
        for out in outgroups:
            if not ali.index(out):
                print '%s not found in the alignment' % out
                return 3
        ali.sort(key=lambda r: 'zzzzzzzz' if r.id in outgroups else r.id)
        AlignmentUtils.save(ali, '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.aln.trimmed.fasta')
        args = dict(plen = (20,40),
                    max_mismatches = 8,
                    min_match_mismatches = 1,
                    first_match_mismatches = 1,
                    first_may_match = 1,
                    AT_first=True,
                    outgroup=len(outgroups))
        fprimers = PrimerFinder.find_discriminating_primers(ali, **args)
        rprimers = PrimerFinder.find_discriminating_primers(ali, reverse=True, **args)
        pairs = PrimerFinder.compile_pairs(fprimers, rprimers, min_prod, 'SSBa')
        if not pairs:
            print '\nNo suitable primer pairs found'
            return 3
        PrimerFinder.print_pairs(pairs)
        orig_ali = PrimerFinder.add_pairs_to_alignment(pairs, orig_ali)
        AlignmentUtils.save(orig_ali, '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.with_primers.aln.fasta')
        print 'Done'


if __name__ == '__main__':
    Main(run=True)
