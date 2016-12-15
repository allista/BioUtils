#!/usr/bin/python
# coding=utf-8
from BioUtils.NCBI import BatchEntrez

from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.SeqUtils import load_files, safe_write, filename_for_record

import re

coded_re = re.compile(r'(complement\()?(\w+\d+\.?\d)\)?')

class Main(MPMain):
    description = 'Fetch source nucleotide sequences that encode provided proteins.'

    def _main(self):
        accessions_file = '/home/allis/Documents/INMI/Geobacillus-COX/CoxL-markers.txt'
        with open(accessions_file) as inp:
            accessions = set(line for line in (l.strip() for l in inp) if line)
        seqs = load_files(self.abort_event,
                          ['/home/allis/Documents/INMI/Geobacillus-COX/CoxL-analysis.files/CoxL-analysis.gb'])
        targets = filter(lambda s: s.id in accessions, seqs)
        to_fetch = set()
        for t in targets:
            for f in t.features:
                coded = f.qualifiers.get('coded_by')
                if not coded: continue
                m = coded_re.match(coded[0])
                if m: to_fetch.add(m.group(2))
        entrez = BatchEntrez(self.abort_event, 'allista@gmail.com')
        recs = entrez.get_records_for_terms(list(to_fetch), 'nucleotide')
        if recs:
            for r in recs:
                fname = filename_for_record(r)
                print 'Saving: %s' % fname
                safe_write(r, '/home/allis/Documents/INMI/Geobacillus-COX/genomes/'+fname)
        else: print 'No records were fetched'
        print 'Done'

if __name__ == '__main__':
    Main(run=True)


