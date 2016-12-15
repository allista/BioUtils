# coding=utf-8

'''
Created on Mar 16, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os

from collections import Counter

from BioUtils.Tools.Misc import ListDB
from BioUtils.Tools.Text import print_table
from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.AlignmentUtils import AlignmentUtils, AlignmentExt
from BioUtils.Tools.Output import user_message


class Main(MPMain):
    class LetterStats(object):
        def __init__(self, s):
            self.letter = ''
            self.num = 0
            self.freq = 0
            self.freq_no_gaps = 0
            #count frequencies
            f = Counter(s.upper())
            lm = ''
            lmn = 0
            for l in 'ATGCU':
                ln = f[l]
                if ln > lmn:
                    lm = l
                    lmn = ln
            if lm:
                clen = len(s)
                self.letter = lm
                self.num = lmn
                self.freq = lmn / float(clen)
                self.freq_no_gaps = lmn / float(clen - f['.'] - f['-'])

        def __str__(self):
            return ('%s %5.1f%% (%5.1f%% without gaps)'
                    % (self.letter, self.freq * 100, self.freq_no_gaps * 100))

    @staticmethod
    def _ref_index(i, ref):
        f = Counter(ref[:i])
        if ref[i] in '.-': print 'Warning, reference has gap at %d position' % (i+1)#test
        return i - f['.'] - f['-']

    @staticmethod
    def _col_index(i, ref):
        _i = 0
        for j in xrange(len(ref)):
            if ref[j] not in '.-': _i += 1
            if _i == i: return j


    def _main(self):
        alifile = '/home/allis/Documents/INMI/Fervidicoccales-signature/arb-silva.de_2016-04-06_id331139.fasta'
        group_names = ['Fervidicoccales', 'Acidilobales', 'Desulfurococcales', 'Thermoproteales:Sulfolobales', 'Other']
        predefined_positions = [34, 501, 544, 1244, 1293]
        ref_name = 'Escherichia'
        reference = None
        groups = ListDB()
        with user_message('Loadding initial alignment...', '\n'):
            ali = AlignmentUtils.load_first(alifile)
            if not ali: return 1
        with user_message('Sorting alignment into subgroups...', '\n'):
            for rec in ali:
                if ref_name in rec.description:
                    reference = rec
                    continue
                found = False
                for g in group_names:
                    for k in g.split(':'):
                        if k in rec.description:
                            groups[g] = rec
                            found = True
                            break
                if not found: groups['Other'] = rec
            groups = dict((n, AlignmentExt(groups[n])) for n in groups)
        ali_len = ali.get_alignment_length()
        predefined_positions = [self._col_index(i, reference) for i in predefined_positions]

        print ('\nReference sequence:\n>%s\n%s' %
               (reference.description, str(reference.seq).replace('.', '').replace('-', '')))
        print '\nAlignment: %d seqs, %d columns' % (len(ali), ali_len)
        print print_table([(g, '%d sequences' % len(groups[g])) for g in group_names])
        print

        main_group = group_names[0]
        main_ali = groups[main_group]
        others = group_names[1:]
        for ci in xrange(ali_len):
            main_letter = self.LetterStats(main_ali[:,ci])
            predef = ci in predefined_positions
            if predef or main_letter.freq_no_gaps >= 0.95 and main_letter.freq > 0.5:
                other_letters = [self.LetterStats(groups[g][:,ci]) for g in others]
                if predef or any(l.letter != main_letter.letter for l in other_letters):
                    print ('------------------ E.coli position: %d ---------------------' %
                           (self._ref_index(ci, reference)+1))
                    print print_table([(main_group, str(main_letter))]+
                                      [(g, str(l)) for g, l in zip(others, other_letters)])
                    print
        print 'Done'

if __name__ == '__main__':
    Main(run=True)
