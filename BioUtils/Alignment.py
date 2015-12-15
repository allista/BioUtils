# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# BioUtils is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# indicator_gddccontrol is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on Jul 20, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import sys
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import MultipleSeqAlignment


class Alignment(object):
    '''
    This class provides abstract interface for sequence alignment using Biopython
    '''

    def __init__(self):
        pass
    #end def
    
    
    def align(self, seq_records):
        '''Align given sequences, return an alignment object'''
        align_cli = MuscleCommandline(clwstrict=True)
        try:
            child = subprocess.Popen(str(align_cli),
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=(sys.platform!="win32"))
            SeqIO.write(seq_records, child.stdin, 'fasta')
            child.stdin.close()
            alignment = AlignIO.read(child.stdout, "clustal")
        except Exception, e:
            raise OSError('Alignment.align: unable to align sequences: ' + e.message)
        return alignment
    #end def
    
    
    def align_translated(self, seq_records, trans_table=1, force_trans_table=False):
        '''Translate given sequences, align them and untranslate.
        Return an alignment object'''
        trans_seqs = self._translate(seq_records, trans_table, force_trans_table)
        alignment  = self.align(trans_seqs)
        untrans_aligned_seqs = self._untranslate(alignment, seq_records)
        return MultipleSeqAlignment(untrans_aligned_seqs)
    #end def
    
    
    def _translate(self, seq_records, trans_table=1, force_trans_table=False):
        '''Translate given sequences using trans_table. If some sequence feature 
        contains a translation table number, use it instead. If force_trans_table is True, 
        always use trans_table.'''
        trans_seqs =  list()
        for record in seq_records:
            #set the alphabet to DNA
            try:
                record.seq.alphabet = IUPAC.ambiguous_dna
            except Exception,e:
                #provided sequences SHOULD be DNA ones
                raise ValueError('Alignment.translate: unable to set alphabet of the sequence:\n%s\n%s' \
                                 % (record.format('fasta'), e.message))        
            #determine translation table
            translation_table = -1
            if force_trans_table:
                #force a translation table
                translation_table = trans_table
            else:
                #see if a translation table is defined in qualifiers
                for feature in record.features:
                    try:
                        translation_table = int(feature.qualifiers['transl_table'][0])
                        break
                    except: pass
                if translation_table < 0: translation_table = trans_table
            #do a translation
            trans_seq = record.seq.translate(table=translation_table, stop_symbol="X")
            trans_seq_rec = SeqRecord(trans_seq, id=record.id)
            trans_seq_rec.name = record.name
            trans_seq_rec.description = record.description
            trans_seqs.append(trans_seq_rec)
        return trans_seqs
    #end def    
    
    
    def _untranslate(self, trans_seq_records, orig_seq_records):
        '''Untranslate given protein sequences using original nucleotide sequences, 
        thus conserving original triplets.'''
        #make a dictionary of sequences
        seq_dict = dict()
        for seq in orig_seq_records:
            seq_dict[seq.id] = seq
        #untranslate sequences
        untrans_seqs = list()
        for record in trans_seq_records:
            orig_seq = str(seq_dict[record.id].seq)
            rev_seq = ""
            l = 0
            for letter in record.seq:
                if letter == "-": rev_seq += "-"*3
                else:
                    if l+3 < len(orig_seq):
                        rev_seq += orig_seq[l:l+3]
                    else: rev_seq += orig_seq[l:]
                    l += 3
            if l < len(orig_seq): rev_seq += orig_seq[l:]
            untrans_seqs.append(SeqRecord(Seq(rev_seq, IUPAC.ambiguous_dna), id=record.id))
            #check results
            if rev_seq.replace("-", "") != orig_seq:
                raise Exception('Alignment.untranslate: resulting sequence differs from the original.')
        return untrans_seqs
    #end def
#end class