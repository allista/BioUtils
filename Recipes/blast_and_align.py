# coding=utf-8

'''
Created on Mar 15, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import re

from BioUtils.Tools.Output import user_message
from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.NCBI import BlastCLI, BlastFilter, BlastID
from BioUtils.SeqUtils import simple_rec
from BioUtils.AlignmentUtils import AlignmentUtils
from BioUtils.PhyloUtils import PhyloUtils
from BioUtils.Taxonomy import Organisms, Lineage

class Main(MPMain):
    def fix_ids(self, seqs):
        for s in seqs:
            desc = s.description.replace(s.id, '').strip()
            s.id = s.name = BlastID.extract(desc) or desc
            s.description = desc
    
    def _main(self):
        query = simple_rec('AAACTGGGGCTAATACCCGATGGGTGAGGAGGCCTGGAATGGTTCTTCACCGAAAAGACGTTGAGACCATGCTTTTCAACGTTGCCTAAGGATGGGGCCGCGTCCGATCAGGTTGTTGGTGGGGTAACGGCTCACCAAGCCTATAACCGGTACGGGCCGTGGGAGCGGAAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTACGGGGCGCAGCAGTCGCGAAAACTCCGCAATGCGCGAAAGCGTGACGGGGCTACCCCGAGTGCCGTCCGCTGAGGATGGCTTTTCCCCGGTGTAATGAGCCTGGGGAATAAGGAGAGGGCAAGCCTGGTGTCAGCCGCCGCGGTAATACCAGCTCTCCGAGTGGTAGGGATGATTATTGGGCTTAAAGCGTCCGTAGCCAGCCCGGCAAGTCTCCCGTTAAATCCAGCGACCTAATCGTTGGGCTGCGGAAGATACTGTTGGGCTAGGGGGCGGGAGAGGCCGACGGTATTCCCGGGGTAGGGGTGAAATCCTATAATCCTGGGAGGACCACCAGTGGCGAAGGCTGTCGGCTAGAACGCGCTCGACGGTGAGGGACGAAAGCTGGGGGAGCGAACTGGATTAGATACCCGGGTAGTCCCAGCTGTAAACGATGCGGGCTAGGTGTTGGGGTGGCTACGAGCCACCTCAGTGCCGCAGGGAAGCCATTAAGCCCGCCGCCTGGGAAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGCGTGAAGCTTGCGGTTTAATTGGAGTCAACGCCGGGAACCTTACCGGGGGCGACAGCAGGATGAGGGCCAGATTGAAGGTCTTGCTTGACAAGCTGAGAGGAGGTGCATGGCCGTCGCCAGTTCGTGCCGTGAGGTGTCCTGTTAAGTCAGGCAACGATCGAGACCCGCACCCTTAGTTGCAACCCCTGCGGAACCCGCAGGGGGCACACTACGGGAACTGCCGCCGATAAGGCGGAGGAAGGAGCGGGCCACGGCAGGTCAGTATGCCCCGAATCCCCCGGGCCACACGCGAGCTGCAATGGCAGAGACAATGGGTTCCAACCTTGAAAGAGGGAGGTAATCCCTAAACCCTGCCTCAGTTGGGATCGAGGGCTGCAACCCGCCCTCGTGAACATGGAATGCCTAGTAATCGCGTGTCATCATCGCGCGGTGAATACGTCCCCGCTCCTTGCACACACCGCCCGTCGCTCCATCCGAGTGGGGTTTGGGTGAGGCGTGGTCTGTTGGCCGCGTCGAATCTAGGCTTCGCGAGGAGGGAGAAGTCGTAACAAGGTGGCCGTAGGGGAACCTGCGGCCGGATCACCTCCT',
                           'BA2-16S')
        suns_db = '/home/allis/Documents/INMI/SunS-metagenome/BlastDB-big/sunspring_meta'
        silva_db = '/home/allis/Documents/INMI/SILVA-DB/SILVA_123_SSURef_Nr99_tax_silva'

        additions = [simple_rec('AAACTGGGGCTAATCCCCCATAGGCCTGGGGTACTGGAAGGTCCCCAGGCCGAAAGGG------GACCGTA-----AGGTCCCGCCCGAGGATGGGCCGGCGGCCGATTAGGTAGTTGGTGGGGTAACGGCCCACCAAG--CCGAAGATCGGTACGGGCC-GTGAGAGCGGGAGCCCGGAGATGGACA---CTGAGACACGGGTCCAGGCCCTACGGGGCGCAGCAGGCGCGAAACC-TCCGCAATGCGGGAAACCGCGACGGGGGGACCCCCAGTGCCGTGCCTCTGGC-----ACGGCTTTTCCGGAGTG-TAAAAAGCTCCGGGAATAAGGGCTGGGCAAGGCCGGTGGC-AGCCGCCGCGGTAATACCGGCGGCCCGAGTGGTGGCCACTATTATTGGGCCTAAAGCGGCCGTAGCCGGGCCCGTAAGTCCCTGGCG-AAATCCCACGGCTCAACCGTGGGGCTCGCTGGGGATACTGCGG-GCCTTGGGACCGGGAGAGGCCGGGGGTACC-CCCGGGGTAGGGGTGAAATCCTATAATCCCGGGGGGACCGCCAGT-GGCGAAGGCGCCC--GGCTGGAACGGGTCCGACGGTGAGGGCCGAAGGCC-AGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCTGGCTGTAAAGGATGCGGGCTAGGTGTCGGGCGAG-CTTCGAGCTCGC-CCGGTGCCGTAGGGAAGCCGTTAAGCCCGCCGCC-TGGGGAGTACGGCCGCAAGGCT-GAAACTTAAAGGAATT-GGCGGGGGAGC-ACTACAAGGGGTGGAGCGTGCGGTTTAATTGGATTCAACGCCGGGAACCTCACCGGGGGCGACGGCAGGATGAA-GGCCAGGCTGAAGGTCTTGCCGGACGCGCCGAGAGGAG-----------------------------------GTGCATGGCCGCCGTCAGCTCGTACCGTGAGGCGTCCA-CTTAAGTGTGGTAACGAGCGAGACCCGC--GCCCCCAGTTGCCAGTCCCTCCCGCTGGGA---GGGAGGC-ACTCTGGGGGG-ACTGCCGGCGAT-AAGCCGGAGGAAGGGGCGGGCGACGGTAGGTCAGTATG-CCCCGAAACCC-CCGGGCT-ACACGCGCGCTACAATGGGCGGGACAATGGGA-CCCGACCCCGAAAGGGGAAGGGAATCCCCTAAACCCGCCCTCAGTTCGGATCGCGGGCTG-CAACTCGCCCGCGTGAAGC-TGGAAT-CCCTAGTACCCGCGCGTCATCATCGCGCGGCGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCACTCCACCCGAG-CGGGGCCC-GGGTGAGGCCCGATCTCCTTCGGGAGGTCGGGTCGAGCCTGGGCTC-CGTGAGGGGGG-AGAAGTCGTAACAAGGTAGCC------------------------------'.replace('-', ''),
                                'Thermococcus_chitonophagus'),
                     simple_rec('AAACTGGGATTAATACCCACTAAATGATAATACCTGGAATGGCTTATCATTGAAAGAC-TCTGGAAACATGCTTC-CAGCGTCGCCCAAGG-------------------------------------------------------------------------------GGAGCCCGGAGATGGAAA---CTGAGACAAGGTTCCAGGCCCTACGGGGCGCAGCAGGCGCGAAACC-TCCACAATGCGCGAAAGCGTGATGGGGTTATCCCGAGTGCCGTCCGATGAGG-----ATGGCTTTTCCTCGGTG-TAAGGATCCGAGGGAATAAAGGGGGGGCAAGACTGGTGTC-AGCCGCCGCGGTAATACCAGCTCCCTGAGTGGTAAGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTATCAAGTCTCTTGTT-AAACCCAGTGATTCAATCATTGACCT-GCAAGAGATACTGTTA-TGCTAGAGGACGGGAGAGGTCGACGG---------GGGTAGGGGTGAAATCCTATAATCCTTGGAGGACCACCAGT-GGCGAAGGCGGTC--GACTAGAACGTGCCTGACGGTGAGGGACGAAAGCT-GGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCCAGCTGTAAACGATGCGGGCTAGGTGTTGGGGTAG-CTACGAGCTACT-CCAGTGCCGCAGAGAAGTTGTTAAGCCCGCCGCC-TGGGGAGTACGGCCGCAAGGCT-GAAACTTAAAGGAATT-GGCGGGGGAGC-ACCACAAGGGGTGAAGGCTGCGGTTTAATTGGAGTCAACGCCGGGAACCTTACCGGGGCTGACAGCAGAGTGAA-GGCCAGACTGAAGATCTTGCCAGACAAGCTGAGAGGAGGTGCATGAAGATCTTGCCAGACAAGCTGAGAGGAGGTGCATGGCCGTCGCCAGTTCGTGCCGTGAGGTGTCCT-GTTAAGTCAGGCAACGAACGAGACCCCC--ACTGTTAGTTGCCAGCGAATTCCAACGGAAT--GTCGGGC-ACACTAACAGG-ACTGCCACCGAT-AAGGTGGAGGAAGGAGGGGGCAACGGCAGGTCAGTATG-CCCC--------------------------------------------------------------------------------------------------------------GAACTCGCCCTCATGAACA-TGGAAT-CCCTAGTAACCGCGTGTCATCATCGCGCGGTGAATACGTCCCCGCTCCTTGCACACACCGCCCGTCGCTCCATCCAAG-TCGGGTCT-AGATGAGGCGCAGTCTTCT-----TGGCTACGTCGAATCTGGGTTC-GGTGAGGGGGG-AGAAGTCGTAACAAGGTGGCCGTAGGGGAACCTGCGGCCGGATCACCTCCT'.replace('-', ''),
                                'SMTZ1-55'),
                     simple_rec('ACTCCGGTTGATCCTGCCGGACCCCACTGCTATCGGGGTAGGACTTAACCATGCGAGTTGTGCGTCCCCAAGCCATGGTGGGGGCGCGGCATACGGCTCAGTAACACGTGGCTAACCTAGCCTTTGGACGGGGACAACCCCGGGAAACTGGGGCTAATCCCCGATGGGTGGGAAGGCCTGGAATGGTTTCCCACCGAAAGGGCGTCTGAACCATGCTTCAGGCGTTGCCGAAGGATGGGGCCGCGGCCGATCAGGTTGTTGGTGAGGTAACGGCTCACCAAGCCTATAACCGGTACGGGCCGTGAGAGCGGGAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTACGGGGCGCAGCAGGTGCGAAAACTCCGCGATGCGCGAAAGCGTGACGGGGCTATCCCGAGTGCCGTCCGCTGAGGATGGCTTTTCCCCGGTGTAGGGAGCCGGGGGAATAAGGAGAGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCTCTCCGAGTGGTGGGGACAATTATTGGGCTTAAAGCGTCCGTAGCCGGCCCATCAAGTCTCTTGTTAAATCCAGCGATCCAATCGCTGGACTGCGGGAGATACTGCTGGGCTAGGGGGCGGGAGAAGCCGATGGTATTCTCGGGGTAGGGGTGAAATCCTATAATCCCGGGAGGACCACCAGTGGCGTAGGCGGTCGGCTAGAACGCGCCCGACGGTGAGGGACGAAAGCTGGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCCAGCCGTAAACGATGCGGGCTAGGTGTTGGGGTGGCTACGAGCCACCCCAGTGCCGCATGGAAGCAATTAAGCCCGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGAAGCTTGCGGTTTAATTGGAGTCAACGCCGGGAAAGGAACAGCGTTTTGTTGTTCCTCTGGATACCTTACCGGGGGCGACAGCAGGATGAAGGCCAGATTGAAGGTCTTGCTGGACGAGCTGAGAGGAGGTGCATGGCCGTCGCCAGTTCGTGCCGTGAGGTGTCCTGTTAAGTCAGGTAACGATCGAGACCCACACCCCCAGTTGCTACCTCTTCGGAGGGCACTCTAGGGGTACTGCCGCCGATAAGGCGGAGGAAGGAGTGGGCCACGGCAGGTCAGTATGCCCCGAATCCCCCGGGCCACACGCGAGCTGCAATGGCAAGGACAATGGGTTCTGACCCCGAGAGGGGAAGGTAATCCCGAAACCCTGCCTCAGTTGGGATCGAGGGCTGAAACCCGCCCTCGTGAACATGGAATCCCTAGTAATCGCGGGTCACCAGCCCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGCTCCATCCGAGTGGGGTTTAGGTGAGGCGTGGTCCTTGTGGCTGTGTCGAATCTAGGCTTCGCGAGGAGGGAGAAGTCGTAACAAGGTGGCCGTAGGGGAACCTGCGGCCGGATCACCTC',
                                'BA1-16S'),
                     simple_rec('CTGGTGGAAATATAGAAGAGGCCAAATCCGGGGTTCAGGCCGCCCGGGGTAATTACCCGTTGTCGGAGTGGGGGGGGGACGCTATTGGGGCTTAAGCCATCGTTAGCCCGTTTGACCAGGTCTCTTGTTAAATCAGGCGGATTTATTGGTCGATTGCAGGAGATTATGTTCGTCTTAGGGGCCGGAGGAGTCAACAGTATTCCCGGGGTAGGAGTGAATGCCTATATTCCCGGAGGTACCACCAGTGGGGACGCCGTTGGTATAGAACGCGCCGGCCGGTGATGGAATGAAAGTGAGGGAACCGACCCGAATTAGATACCGGGGTATTGCTACCGTTAACCGATGCAGCTTAGGTGTTCGGGTGGTTACTAGCCATTCGAGTGCGCCAGGGAAGCTGTCAGGCTTACCGCTTGGGAAGTGCGGCTGCAGGGCCAAAACTTAAGGAAATCGCCGGGGAAGCACCCCAGGGGGTGAAGCTTGCGCTTTAATGGAATTCACCGCGGTAATTCTCACCGGGGGAGCCACCAGGAGGAAAGCCAGATTAAAGTTCTTGTTGGCGGAGTGGAGAGGAGGTGCATGCCGTTCGCCAGTTCTTCCCGGGAGGTTCTTGTTAGTTCAGCCACCGATGAGGACCGCCATCCCCTGTTGTTATTGGCCTTGCGCCAGGCACACTGGGGAGACCGCCGCCGATAAGGCGGAGGAAGGAGCGGGCCACGGCAGGTCAGTATGCCCCGAATCCCCCGTCCACACGCGAGGGGCAATG',
                                '155a'),
                     simple_rec('CAAGTCCTATAACCGGTACGGGCCGTGGGAGCGGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTACGGGGCGCAGCAGTCGCGAAACCTCCGCAATACGCGAAAGCGTGACGGGGTCATCCCGAGTGCCGTCCGCTGAGGATGGCTTTTCCCCAGTGTAGACAGCTGGGGGAATAAGGAGAGGGCAAGTCGGGTGTCAGCCGCCGCGGTAATACCCGCTCTCCGAGTGGTGGGGACGCTTATTGGGCCTAAAGCATCCGTAGCCGGCTGGACAAGTCCCCTGTTAAATCCAGCGATTTAATCGTTGGACTGCGGGGGATACTGTCCGGCTAGGGGGCGGGAGAGGCCGACGGTATTTCCGGGGTAGGGGTGAAATCCTATAATCCCGGGAGGACCACCAGTGGCGAAGGCTGTCGGCTAGAACGCGCCCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCCAGCCGTAAACGATGCAGGCTAGGTGTTCGGGTGGCTACGTGCCACTCGAGTGCCGCAGGGAAGCTGTTAAGCCTGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGAAGCTTGCGGTTTAATTGGAGTCAACGCCGGAAATCTCACCGGGGGAGACAGCAGGATGAAAGCCAGATTAAAGGTCTTGCTAGACGAGCTGAGAGGAGGTGCATGGCCGTCGCCAGTTCGTGCCGTGAGGTGTCCTGTTAAGTCAAGGCAACGATCGAGACTCGCATCCTCTGTTGCTACTACCCTTGCGCCAAGGCACACTGGGGGAGACCGCCGCTCGATAAGGCGGAAGGAAGGAGCGGCCCACGGCAGTCAGTATGCCCCGAATTCCCTCGGCCACACGCAAGCTGCAATG',       
                                '156a'),
                     simple_rec('GGGGATCGGGGCATACTGACCTGCCGTGGCCCGCTCCTTCCTCCGCCTTATCGGCGGCGGTCTCCCCAGTGTGCCTGGCGCAAGGGCAGTAACAACAGGGGATGGGGGTCTCGATCGGTGGCTGGCTTAACAGGAAACCTCACGGGACGAACTGGCGAACGGCATGGACCTTCTCTCAACTTGGCTAAGAAGAACTTTAATCTGGCTTTCATTCTGGTGGCTTCCCCGGTGAGAATTCCGGCGGTGACTCCCAATAAAACGCAAGCTTCACCCCTTGGGGTGGTTCCCCGGCCATTTCTTTAAGGTTCAAGCTTTGCGGCGGTATTCCCAAGCGGCAAGGTTAACAGCTTCCCTGCCGCACTCGAGTGGCACGTAACCACCCGAACAACTAACCTGCATCCGTTACCGGTTGGACTAACCCGGTATCTAATCCGGGTCGCTCCCCCAGCCTTCATTCCTTCACCGTCCGGCGCGGTTCTAAGCGACCGGCTTTCGCACTTGTGGTTCCTCCCGGGGATTATAAGAATTCACCCCTACCCCGGAAATTACGGTCCGGCTCCTCCGGCCCCTAACCCGACACGTAATCCCCCGCCAGTTCAACCGATTAAATCCGCTTGAATTTAACAAGGGGGACCTTGTCCAGCCGGCCTACGGATGCTTTAAGGCCCAATAAGCCGTCCCCACCACTCCGAGAGCGGGTAATAACCGCGGCCGGCCTGACAACCGACCTGGCCTCTCCTAAATCCCCCAGCTGTTCACACTTGGGAAAGGGCATTCCTCAGCGAACGGCACTTCGGGATGAACCCGTCACGCTTTCGCGTAATTGCGGGAAGGTTTCGCGAACTGCTGCGCCCCGTAAAGGCCTGGGTCCTTGTGTCTCAAATTGCCCCATCTCCGGGCTATACGCTCTCCACGGGCCCGTACC', 
                                '157a')
                     ]
        #prepare filter
        filt = BlastFilter(lambda a, r: a.hsps[0].align_length > 1100)
        filt.AND = BlastFilter(lambda a, r: all(hsp.score > 500 for hsp in a.hsps))
        filt.AND.AND = BlastFilter(lambda a, r: all(hsp.identities/float(hsp.align_length) > 0.8 for hsp in a.hsps))
        #make ring-blast
        blast = BlastCLI(self.abort_event)
        orig_seqs = blast.ring_blast(query, suns_db, 100, filt, 3)
        if not orig_seqs:
            print 'No blast results.'
            return 1
        nseqs = len(orig_seqs)
        print 'RingBlast to:\n%s\nreturned %d sequences.\n' % (suns_db, nseqs)
        #save an initial alignment
        self.fix_ids(orig_seqs)
        alifile = '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.aln.fasta'
        with user_message('Aligning retrieved sequences...', '\n'):
            if not AlignmentUtils.align(orig_seqs+[query]+additions, outfile=alifile): return 3
        #search for additional homologs
        add_seqs = blast.ring_blast(orig_seqs, silva_db, 100, filt, 0)
        if add_seqs:
            self.fix_ids(add_seqs)
            print 'RingBlast to:\n%s\nreturned %d additional sequences.\n' % (silva_db, len(add_seqs))
        #build an alignment
        seqs = orig_seqs+add_seqs+[query]+additions
        alifile = '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.big.aln.fasta'
        with user_message('Aligning retrieved sequences...', '\n'):
            if not AlignmentUtils.align(seqs, outfile=alifile): return 3
        #build a tree 
        treefile = '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.big.aln.tre'
        if not PhyloUtils.build_fast_tree(alifile, treefile): return 4
        #annotate the tree
        if False:
            with open('/home/allis/Documents/INMI/16S/SSBaF4-SSBaR4-1_243072232-iPCR-report.txt') as inp:
    #            SSBaF4-SSBaR4_65397396-iPCR-report.txt
                sids = set()
                len_re = re.compile(r'(\s|^)(\d+)(\sbp|\\s*:)?', re.MULTILINE)
                entry = False
                cur_sid = None
                cur_len = -1
                for l in inp:
                    if l == '========= histograms and electrophorograms of PCR products of each hit =========': break
                    if l.startswith('---'): 
                        entry = False
                        if cur_sid and cur_len > 0 and abs(cur_len-920) < 60:
                            sids.add(cur_sid)
                        cur_sid = None
                        cur_len = -1
                        continue
                    if entry or '#' in l:
                        entry = True
                        plen = len_re.search(l)
                        if plen: cur_len = int(plen.group(2))
                        sid = BlastID.extract(l)[0]
                        if sid: cur_sid = sid
        organisms = Organisms.from_records(seqs)
        if PhyloUtils.annotate_tree(treefile, organisms, 
                                    reroot_at='Thermococcus_chitonophagus',
#                                    beautify_leafs=True,
#                                    collapse_taxa=['miscellaneous crenarchaeotic group', 'thaumarchaeota'],
#                                    collapse_last=True,
#                                    collapse_hard=True,
#                                    mark_leafs=sids,
                                    mark_leafs=[r.id for r in orig_seqs+[query]+additions],
                                    lineage_colors={'miscellaneous crenarchaeotic group':(0, 0, 255),
                                                    'thaumarchaeta':(255,0,0)},
                                    top_lineage=Lineage('archaea')): return 0
        return 2

if __name__ == '__main__':
    Main(run=True)
