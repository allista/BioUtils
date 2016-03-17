# coding=utf-8

'''
Created on Mar 16, 2016

@author: Allis Tauri <allista@gmail.com>
'''

from BioUtils.Tools.Text import issingleletter
from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.AlignmentUtils import AlignmentUtils, AlignmentExt
from BioUtils.SeqUtils import copy_attrs

from DegenPrimer.Primer import Primer

class Main(MPMain):
    def _match_mismatch(self, *s):
        if not s: return False
        ret = True
        for si in s:
            ret &= issingleletter(si[:-1]) and not issingleletter(si)
            if not ret: break
        return ret

    def _find_primers(self, ali,
                      plen = 25,
                      max_mismatches = 4,
                      min_match_mismatches = 3,
                      first_matches = 2):
        primers = []
        for s in list(ali.cols)[:-plen]:
            mis = 0
            mmis = min_match_mismatches
            pali = ali[:,s:s+plen]
            if pali[1,-1].upper() not in ('A', 'T'):
                continue
            good = True
            for i in pali.rcols:
                col = pali[:,i]
                if plen-i <= first_matches:
                    if not self._match_mismatch(col):
                        good = False
                        break
                    continue
                if '-' in col[:-1]:
                    good = False
                    break
                if self._match_mismatch(col): mmis -= 1
                if not issingleletter(col[:-1]): mis += 1
            good &= mmis <= 0 and mis <= max_mismatches
            if not good: continue
            print 'position: %d\n%s\n' % (s, str(pali))
            primers.append((s+plen, pali))
        return primers

    def _main(self):
        min_prod = 400
        outgroup_name = 'gnl|BL_ORD_ID|72134:c1348-71'
        exclude = ['NR_042737', 'SMTZ1-55', 'BA2-16S']
        alifile = '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.aln.fasta'
        alis =  AlignmentUtils.load(alifile)
        if not alis: return 1
        orig_ali = AlignmentExt.from_msa(alis[0])
        ali = orig_ali.remove(*exclude).trim()
        outgroup_index = ali.index(outgroup_name)
        if outgroup_index is None:
            print 'Outgroup is not found in the alignment'
            return 2
        ali.sort(key=lambda r: r.id if r.id != outgroup_name else 'zzzzzzzz')
        ali_len = ali.get_alignment_length()
        outgroup_index = ali_len-1
        AlignmentUtils.save(ali, '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.aln.trimmed.aln')
        fprimers = self._find_primers(ali)
        rprimers = self._find_primers(ali.reverse_complement())
        pairs = []
        for i, (fs, fp) in enumerate(fprimers):
            start = fs
            fprimer = Primer.from_sequences(fp[:-1], 1, 'SSBaF%d' % (i+1))
            for j, (rs, rp) in enumerate(rprimers):
                end = ali_len-rs
                print start, end, min_prod
                if end-start <= min_prod: continue
                pairs.append((fprimer, Primer.from_sequences(rp[:-1], 1, 'SSBaR%d' % (j+1))))
        if not pairs:
            print 'No suitable primer pairs found'
            return 3
        added = set()
        for i, (fp, rp) in enumerate(pairs):
            print 'pair %d' % (i+1)
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
        orig_ali = AlignmentUtils.align(orig_ali)
        AlignmentUtils.save(orig_ali, '/home/allis/Documents/INMI/SunS-metagenome/Bathy/BA2_SunS_16S.with_primers.aln')
        print 'Done'


if __name__ == '__main__':
    Main(run=True)
