# coding=utf-8

'''
Created on Mar 16, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os

from BioUtils.Tools.Text import issingleletter, FilenameParser
from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.AlignmentUtils import AlignmentUtils
from BioUtils.SeqUtils import copy_attrs, SeqView
from BioUtils.Tools.Output import user_message

from DegenPrimer.Primer import Primer

class Main(MPMain):
    def _match_mismatch(self, s, outgroup=1):
        if not s: return False
        return (issingleletter(s[:-outgroup]) 
                and not s[0] in s[-outgroup:]) 

    def _find_primers(self, alignment, 
                      plen = 25,
                      max_mismatches = 3,
                      min_match_mismatches = 3,
                      first_match_mismatches = 3,
                      first_may_match = 0,
                      AT_first = True,
                      outgroup=1):
        primers = []
        if isinstance(plen, (tuple, list)):
            min_len = plen[0]
            max_len = plen[1]
        else: min_len = max_len = plen
        for s in alignment.cols:
            if s < min_len-1: continue
            mlen = min_len
            mis  = 0
            mmis = min_match_mismatches
            pali = alignment[:,max(0,s-max_len):s]
            if AT_first and pali[0,-1].upper() not in 'AT': continue
            good = True
            cur_len = 0
            first_matched = 0
            for i in pali.rcols:
                cur_len += 1
                col = pali[:,i]
                group = col[:-outgroup]
                if cur_len <= first_match_mismatches+first_matched:
                    if cur_len <= first_may_match:
                        if not issingleletter(group):
                            good = False
                            break
                        elif not self._match_mismatch(col, outgroup):
                            first_matched += 1
                    elif not self._match_mismatch(col, outgroup):
                        good = False
                        break
                    continue
                if '-' in group:
                    if issingleletter(group, '-'): mlen += 1
                    else:
                        good = False
                        break
                if self._match_mismatch(col, outgroup): mmis -= 1
                if not issingleletter(group): mis += 1
                if cur_len >= mlen and mmis <= 0 and mis <= max_mismatches: break
            good &= cur_len >= mlen and mmis <= 0 and mis <= max_mismatches
            if not good: continue
            primers.append((s+1, pali[:-outgroup,pali.get_alignment_length()-cur_len:].trim()))
            print '\nposition: %d\n%s' % primers[-1]#test
        return primers

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
        ali_len = ali.get_alignment_length()
        AlignmentUtils.save(ali, '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.aln.trimmed.fasta')
        args = dict(plen = (20,40),
                    max_mismatches = 8,
                    min_match_mismatches = 1,
                    first_match_mismatches = 1,
                    first_may_match = 1,
                    AT_first=True,
                    outgroup=len(outgroups))
        fprimers = self._find_primers(ali, **args)
        rprimers = self._find_primers(ali.reverse_complement(), **args)
        pairs = []
        for i, (fs, fp) in enumerate(fprimers):
            start = fs
            fprimer = Primer.from_sequences(fp[:-1], 1, 'SSBaF%d' % fs)
            for _j, (rs, rp) in enumerate(rprimers):
                end = ali_len-rs
                if end-start <= min_prod: continue
                pairs.append((fprimer, Primer.from_sequences(rp[:-1], 1, 'SSBaR%d' % (ali_len-rs+1))))
        if not pairs:
            print '\nNo suitable primer pairs found'
            return 3
        added = set()
        for i, (fp, rp) in enumerate(pairs):
            print '\npair %d' % (i+1)
            print '%s: %s' % (fp.id, fp)
            print '%s: %s' % (rp.id, rp)
            if fp.id not in added:
                orig_ali.append(fp.master_sequence+'-'*(orig_ali.get_alignment_length()-len(fp)))
                added.add(fp.id)
            if rp.id not in added:
                orig_ali.append(copy_attrs(rp.master_sequence,
                                           rp.master_sequence.reverse_complement())+
                                '-'*(orig_ali.get_alignment_length()-len(rp)))
                added.add(rp.id)
        print
        orig_ali = AlignmentUtils.align(orig_ali)
        AlignmentUtils.save(orig_ali, '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.with_primers.aln.fasta')
        print 'Done'


if __name__ == '__main__':
    Main(run=True)
