# coding=utf-8

'''
Created on Mar 16, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import sys
from BioUtils.AlignmentUtils import AlignmentUtils

from BioUtils.SeqUtils import copy_attrs
from BioUtils.Tools.Text import issingleletter
from DegenPrimer.Primer import Primer


class PrimerPair(object):
    def __init__(self, fprimer, rprimer):
        self.forward = fprimer
        self.reverse = rprimer

    def __str__(self):
        return '%s: %s\n%s: %s' % (self.forward.id, self.forward,
                                   self.reverse.id, self.reverse)

class PrimerFinder(object):
    """A set of routines for searching possible primer paris in a nucleotide alignment."""

    @staticmethod
    def _match_mismatch(s, outgroup=1):
        if not s: return False
        return (issingleletter(s[:-outgroup])
                and not s[0] in s[-outgroup:])

    @staticmethod
    def find_specific_primers(alignment,
                              plen=25,
                              max_mismatches=3,
                              min_first_matches=3,
                              AT_first=True,
                              reverse=False):
        primers = []
        if reverse: alignment = alignment.reverse_complement()
        ali_len = alignment.get_alignment_length()
        if isinstance(plen, (tuple, list)):
            min_len = plen[0]
            max_len = plen[1]
        else: min_len = max_len = plen
        for s in alignment.cols:
            if s < min_len - 1: continue
            mlen = min_len
            mis = 0
            pali = alignment[:, max(0, s - max_len):s]
            if AT_first and pali[0, -1].upper() not in 'AT': continue
            good = True
            cur_len = 0
            for i in pali.rcols:
                cur_len += 1
                col = pali[:, i]
                if cur_len <= min_first_matches:
                    if not issingleletter(col):
                        good = False
                        break
                    continue
                if '-' in col:
                    if issingleletter(col, '-'): mlen += 1
                    else:
                        good = False
                        break
                if not issingleletter(col): mis += 1
                if mis > max_mismatches:
                    if cur_len > mlen: cur_len -= 1
                    else: good = False
                    break
            if not good: continue
            pali_len = pali.get_alignment_length()
            if not issingleletter(pali[:, pali_len - cur_len]): cur_len -= 1
            primers.append((ali_len-s if reverse else s+1,
                            pali[:, pali_len - cur_len:].trim()))
        return primers

    @classmethod
    def find_discriminating_primers(cls, alignment,
                                    plen=25,
                                    max_mismatches=3,
                                    min_match_mismatches=3,
                                    first_match_mismatches=3,
                                    first_may_match=0,
                                    AT_first=True,
                                    outgroup=1,
                                    reverse=False):
        primers = []
        if reverse: alignment = alignment.reverse_complement()
        ali_len = alignment.get_alignment_length()
        if isinstance(plen, (tuple, list)):
            min_len = plen[0]
            max_len = plen[1]
        else: min_len = max_len = plen
        for s in alignment.cols:
            if s < min_len - 1: continue
            mlen = min_len
            mis = 0
            mmis = min_match_mismatches
            pali = alignment[:, max(0, s - max_len):s]
            if AT_first and pali[0, -1].upper() not in 'AT': continue
            good = True
            cur_len = 0
            first_matched = 0
            for i in pali.rcols:
                cur_len += 1
                col = pali[:, i]
                group = col[:-outgroup]
                if cur_len <= first_match_mismatches + first_matched:
                    if cur_len <= first_may_match:
                        if not issingleletter(group):
                            good = False
                            break
                        elif not cls._match_mismatch(col, outgroup):
                            first_matched += 1
                    elif not cls._match_mismatch(col, outgroup):
                        good = False
                        break
                    continue
                if '-' in group:
                    if issingleletter(group, '-'): mlen += 1
                    else:
                        good = False
                        break
                if cls._match_mismatch(col, outgroup): mmis -= 1
                if not issingleletter(group): mis += 1
                if cur_len >= mlen and mmis <= 0 and mis <= max_mismatches: break
            good &= cur_len >= mlen and mmis <= 0 and mis <= max_mismatches
            if not good: continue
            pos = ali_len - s if reverse else s + 1
            primers.append((pos, pali[:-outgroup, pali.get_alignment_length() - cur_len:].trim()))
        return primers

    @staticmethod
    def compile_primers(primer_alis, name_base, reverse=False):
        primers = []
        if reverse:
            for pos, ali in primer_alis:
                primers.append(Primer.from_sequences(ali, 1, '%sR%d' % (name_base, pos)))
        else:
            for pos, ali in primer_alis:
                primers.append(Primer.from_sequences(ali, 1, '%sF%d' % (name_base, pos)))
        return primers

    @staticmethod
    def compile_pairs(falis, ralis, min_len, name_base):
        pairs = []
        for _i, (fs, fp) in enumerate(falis):
            fprimer = Primer.from_sequences(fp, 1, '%sF%d' % (name_base, fs))
            for _j, (rs, rp) in enumerate(ralis):
                if rs - fs <= min_len: continue
                pairs.append(PrimerPair(fprimer, Primer.from_sequences(rp, 1, '%sR%d' % (name_base, rs))))
        return pairs

    @classmethod
    def all_pairs(cls, falis, ralis, name_base):
        return cls.compile_pairs(falis, ralis, -sys.maxint, name_base)

    @staticmethod
    def print_pairs(pairs):
        for i, p in enumerate(pairs):
            print '\npair %d\n%s' % (i + 1, p)

    @staticmethod
    def add_primers_to_alignment(primers, alignment, reverse=False):
        alignment = alignment.clone()
        ali_len = alignment.get_alignment_length()
        if reverse:
            for p in primers:
                alignment.append(copy_attrs(p.master_sequence,
                                            p.master_sequence.reverse_complement()) +
                                 '-' * (ali_len - len(p)))
        else:
            for p in primers:
                alignment.append(p.master_sequence + '-' * (ali_len - len(p)))
        return AlignmentUtils.align(alignment)

    @staticmethod
    def add_pairs_to_alignment(pairs, alignment):
        added = set()
        alignment = alignment.clone()
        ali_len = alignment.get_alignment_length()
        for i, p in enumerate(pairs):
            if p.forward.id not in added:
                alignment.append(p.forward.master_sequence + '-' * (ali_len - len(p.forward)))
                added.add(p.forward.id)
            if p.reverse.id not in added:
                alignment.append(copy_attrs(p.reverse.master_sequence,
                                            p.reverse.master_sequence.reverse_complement()) +
                                 '-' * (ali_len - len(p.reverse)))
                added.add(p.reverse.id)
        return AlignmentUtils.align(alignment)
