
'''
Created on Dec 10, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import choice


def generate_diff_sequences(seq_len, diff_by, seqs_number):
    '''Generate seqs_number sequences of seq_len length which deffer 
    at least by diff_by letters'''
#        maximum_different  = max_combs(diff_by, seq_len)
#        if seqs_number > maximum_different:
#            print '\nWarning: requested number of sequences to generate ' \
#                  'exceeds maximum possible number: %d' % maximum_different
#            #seqs_number = maximum_different
    sequences = []
    bad_seqs  = 0
    max_bad   = 4**(5+(seq_len-diff_by))
    while len(sequences) < seqs_number:
        #generate new sequence
        new_seq = ''.join(choice('ATGC') for i in xrange(seq_len))
        #check differences
        differ = True
        for seq in sequences:
            diffs = 0
            for i in xrange(seq_len):
                diffs += int(seq[i] != new_seq[i])
                if diffs >= diff_by: break
            if diffs < diff_by:
                differ = False 
                break
        #add new sequence if it deffers from each of previously added
        if differ: 
            sequences.append(new_seq)
            bad_seqs = 0
        else: bad_seqs +=1
        if bad_seqs > max_bad: break
    return sequences
#end def


def generate_tags(tag_len, tag_diff_by, num_tags, sequences):
    tags = generate_diff_sequences(tag_len, tag_diff_by, num_tags)
    tagged_sequences = []
    for seq_rec in sequences:
        tagged_sequences.append([SeqRecord(Seq(tag+str(seq_rec.seq), seq_rec.seq.alphabet), 
                                           id = seq_rec.id+('_%d' % tag_index),
                                           name = seq_rec.name+('_%d' % tag_index),
                                           description = seq_rec.description
                                           ) for tag_index, tag in enumerate(tags)])
    return tagged_sequences
#end def
    