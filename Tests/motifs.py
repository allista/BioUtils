#!/usr/bin/python
# coding=utf-8

'''
Created on Jan 6, 2016

@author: Allis Tauri <allista@gmail.com>
'''

from Bio.NeuralNetwork.Gene.Schema import SchemaFinder
from BioUtils.SeqUtils import load_files, load_dir, simple_feature
from BioUtils.Tools import MPMain
import sys, csv

class Main(MPMain):
    genomedir = u'/home/allis/Dropbox/Science/Микра/Thermococcus/sequence/GenBank/Thermococcus/'
    
    promoters = {'Thermococcus_barophilus_Ch5':simple_feature(1500426, 1500139-1),
           'CP000855.1':simple_feature(1434193-1, 1434542), #NA1
           'Thermococcus_barophilus_DT4':simple_feature(1357510, 1357226-1),
           'CP006965.1':simple_feature(1404107, 1403821-1), #ES1
           'CP001398.1':simple_feature(60168, 60417-1), #Tgam
           }
    
    def _main(self):
        genomes = load_dir(self.abort_event, self.genomedir, 'gb', '.*\.gb$')
        if not genomes: return 1
        
        pseqs = {}
        
        for g in genomes:
            if g.id not in self.promoters: continue
            if g.id in pseqs: continue
            pseqs[g.id] = self.promoters[g.id].extract(g)
            print g.id
            print pseqs[g.id]
            print 

        finder = SchemaFinder()
        repo = finder.find(pseqs.values())
        with open('FC-promoter.mfs', 'w') as out:
            writer = csv.writer(out)
            writer.writerow(('motif', 'count')) 
            for pi in repo.get_all():
                writer.writerow((pi, repo.count(pi)))

if __name__ == '__main__':
    main = Main()
    sys.exit(main())