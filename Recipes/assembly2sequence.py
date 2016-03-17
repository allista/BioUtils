# coding=utf-8

'''
Created on Mar 16, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os
from Bio.Alphabet import generic_dna, generic_protein
from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.SeqUtils import simple_rec, load_dir, pretty_rec_name, safe_write, copy_attrs

class Main(MPMain):
    def _main(self):
#        self.argument('path', metavar='path',
#                      type=str, nargs='+',
#                      help='Path to a file with the assembly or a directory of such files.')
#        self.argument('--datatype', metavar='str', default='',
#                            type=str, nargs=1, choices=GapEncoder.datatypes.keys(),
#                            help='Matrix datatype. May be %s.' % ', '.join(GapEncoder.datatypes.keys()))
#        self.argument('--schema', metavar='str',
#                            type=str, nargs=1, choices=GapEncoder.schemas, default=[None],
#                            help='Alignment file schema. May be %s.' % ', '.join(GapEncoder.schemas))
#        self.parse_args()
        
        os.chdir('/home/allis/Documents/INMI/SunS-metagenome/Bathy/genomes')
        allseqs = load_dir(self.abort_event, '.',
                           namefilter='.*\.gbff$', flatten=False)
        if not allseqs: 
            print 'No files were loaded'
            return 1
        for seqs in allseqs:
            if len(seqs) < 2:
                print '%s is represented by a single sequence. Skipping.' % pretty_rec_name(seqs[0]) 
                continue
            seq = simple_rec('', '_id', alphabet=generic_dna)+seqs[0]
            for s in seqs[1:]:
                seq += simple_rec('N'*20, 'gap', feature='gap', alphabet=generic_dna)
                seq += s
            copy_attrs(seqs[0], seq)
            filename = seq.name+'.merged.gb'
            safe_write(seq, filename)
            print '%s saved.' % filename
        print 'Done'
            

if __name__ == '__main__':
    Main(run=True)