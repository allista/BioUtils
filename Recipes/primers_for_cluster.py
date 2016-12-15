#!/usr/bin/python
# coding=utf-8

"""
Created on Mar 19, 2016

@author: Allis Tauri <allista@gmail.com>
"""

from __future__ import print_function
import os

from Bio.Blast.NCBIXML import parse

from BioUtils.NCBI import BlastCLI, BlastFilter, BlastWWW
from BioUtils.SeqUtils import SeqLoader, SeqView, get_indexes_of_genes, safe_write, pretty_rec_name, \
    get_feature_indexes_by_qualifier
from BioUtils.AlignmentUtils import AlignmentUtils
from BioUtils.PrimerFinder import PrimerFinder
from BioUtils.Tools.Multiprocessing import MPMain


def find_primers(segments, foi, pfinder_args, reverse=False, qualifier='ugene_name'):
    genes = []
    for s in segments:
        i = get_feature_indexes_by_qualifier(s, qualifier, foi)
        if not i: continue
        genes.append(s.features[i[0]].extract(s))
    ali = AlignmentUtils.align(genes)
    if not ali: return None
    primer_alis = PrimerFinder.find_specific_primers(ali, reverse=reverse, **pfinder_args)
    return PrimerFinder.compile_primers(primer_alis, foi+'_', reverse), ali


class Main(MPMain):
    def _main(self):
        email = 'allista@gmail.com'
        genome_dir = '/home/allis/Dropbox/Science/Микра/Thermococcus/sequence/GenBank/Thermococcales/Thermococcus/'
        genome = 'Thermococcus_barophilus_Ch5.gb'
        gene = 'TBCH5v1_1369' #cooS
        database = 'nr'
        segment = [3200, 12000]

        seq = SeqLoader.load_file(os.path.join(genome_dir, genome))
        if not seq: raise RuntimeError('No genome loaded')
        seq = seq[0]

        index = get_indexes_of_genes(seq, gene)
        if not index: raise RuntimeError('No gene found')

        feature = seq.features[index[0]]
        query = feature.extract(seq)

        segments_file = 'CO-clusters.gb'
        #get cluster variants if needed
        if not os.path.isfile(segments_file):
            blast_file = 'blast.results.xml'
            if os.path.isfile(blast_file): blast = list(parse(open(blast_file)))
            else: blast = BlastCLI.blast_seq(query, database, 100, remote=True, task='blastn',
                                             parse_results=True, save_results_to='blast.results.xml')
            if not blast: raise RuntimeError('Blast returned no results')
            flt = BlastFilter(lambda hsp, r: hsp.align_length > 700, filter_hsps=True)
            flt(blast)
            queries = []
            for ali in BlastCLI.iter_alignments(blast):
                q = BlastCLI.Query(ali, 'hsp', start_offset=segment[0], end_offset=segment[1])
                if q: queries.append(q)
                print(queries[-1])

            segments = BlastWWW.fetch_queries(email, queries)
            safe_write(segments, segments_file)
            for r in segments: print('[%s] %s: %dbp' % (r.id, pretty_rec_name(r), len(r)))
            return 0

        #find primers in alignments of the selected features
        local_files = [os.path.join(genome_dir, f) for f in
                       ('Thermococcus_barophilus_DT4-complete-genome.gb',
                        'Thermococcus_ST-423.gb',
                        'Thermococcus_CH1-complete.gb')]
        loader = SeqLoader(self.abort_event)
        segments = loader.load_files([segments_file]+local_files)
        fprimers, transF_ali = find_primers(segments, 'transF',
                                            dict(plen=(20, 30),
                                                 max_mismatches=5,
                                                 min_first_matches=3,
                                                 AT_first=True))
        rprimers, cooS_ali = find_primers(segments, 'cooS',
                                          dict(plen=(20, 30),
                                               max_mismatches=4,
                                               min_first_matches=3,
                                               AT_first=True),
                                          reverse=True)
        if not fprimers:
            print('\nNo forward primers found')
            return 1
        if not rprimers:
            print('\nNo reverse primers found')
            return 1
        print('\nForward primers:')
        for p in fprimers: print('%s: %s' % (p.id, p))
        print('\nReverse primers:')
        for p in rprimers: print('%s: %s' % (p.id, p))
        print()
        #add primers to alignments and save them
        transF_ali = PrimerFinder.add_primers_to_alignment(fprimers, transF_ali)
        cooS_ali = PrimerFinder.add_primers_to_alignment(rprimers, cooS_ali, reverse=True)
        AlignmentUtils.save(transF_ali, 'transF.aln')
        AlignmentUtils.save(cooS_ali, 'cooS.aln')

if __name__ == '__main__':
    Main(run=True)
