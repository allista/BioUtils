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
Created on Aug 1, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from copy import deepcopy
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from BioUtils.CommonTools import print_exception

class BlastFetcher(object):
    '''Make a blast search possibly restricted by Entrez query, then fetch
    results as sequences with annotations'''

    def __init__(self, job_id):
        self._job_id = job_id
        self._sequences_filename = self._job_id + '-blast-seqs.gb'
        self._results_filename   = self._job_id + '-blast.xml'
        self._blast_results = None
        self._blast_results_sequences = None
        self.have_results   = False
    #end def
    
    
    def blast(self, query, expect=10, megablast=False, entrez_query=''):
        if not query: return False
        print '\nStarting a BLAST search. This may take awhile...'
        try:
            blast_results = NCBIWWW.qblast('blastn', 
                                           self.database, 
                                           self._query.format('fasta'),
                                           expect       = expect, 
                                           megablast    = megablast, 
                                           entrez_query = entrez_query)
            #save results to a file
            results_file = open(self._results_filename, 'w')
            results_file.write(blast_results.read())
            results_file.close()
            blast_results.close()
            print '\nBlast output was written to:\n   ' + self._results_filename
            #parse results
            results_file  = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            #clear fetched sequences
            self._blast_results_sequences = None
        except Exception, e:
            print 'BlastFetcher.blast: failed to obtain BLAST search results from NCBI.'
            print_exception(e)
            return False
        self.have_results = True
    #end def
    
    
    def load_results(self):
        #load blast results
        print '\nLoading previously saved BLAST results:\n   ', self._results_filename
        try:
            results_file        = open(self._results_filename, 'r')
            self._blast_results = list(NCBIXML.parse(results_file))
            results_file.close()
            self.have_results = True
        except Exception, e:
            print '\nFailed to load blast results.\n', e._message
    #end def
    
    
    def fetch_results(self, email):
        if not self.have_results: return
        #if sequences were already fetched, just return them
        if self._blast_results_sequences:
            return deepcopy(self._blast_results_sequences)
        if not email:
            raise Exception('\nYou should always provide a valid e-mail '
                            'to NCBI when performing an Entrez search.')
        print '\nFetching records from Entrez database. This may take awhile...'
        Entrez.email = email
        #search the ID of each blast result, then fetch corresponding part of the sequence
        self._blast_results_sequences = list()
        for record in self._blast_results:
            for alignment in record.alignments:
                SQID = alignment.hit_id
                try:
                    entrez_ID = Entrez.esearch(db='nucleotide', term=('%s[SQID]' % SQID))
                except Exception, e:
                    print '\nUnable to get Entrez ID for sequence %s - %s' % (SQID, alignment.hit_def)
                    print 'Exception message is:', e.message
                    continue  
                for hsp in alignment.hsps:
                    seq_start = hsp.sbjct_start
                    seq_end   = hsp.sbjct_end
                    try:
                        ehandle   = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', 
                                                  seq_start=str(seq_start), seq_end=str(seq_end),
                                                  id=entrez_ID)
                        seq_record = SeqIO.read(ehandle, format='gb', alphabet=IUPAC.ambiguous_dna).next()
                    except Exception, e:
                        print '\nUnable to fetch sequence: %s - %s [%d:%d]' \
                        % (SQID, alignment.hit_def, seq_start, seq_end)
                        print 'Exception message is:', e.message
                        continue
                    self._blast_results_sequences.append(seq_record)
        #write fetched sequences to a file
        if self._blast_results_sequences:
            try:
                sequences_file      = open(self._sequences_filename, 'w')
                SeqIO.write(self._blast_results_sequences, sequences_file, 'gb')
                sequences_file.close()
            except Exception, e:
                print '\nFailed to save fetched sequences.'
                print  'Exception message is:', e.message
        return deepcopy(self._blast_results_sequences)
    #end def
    
    
    def load_sequences(self):
        #load fetched sequences
        print '\nLoading previously saved sequences:\n   ', self._sequences_filename
        try:
            sequences_file      = open(self._sequences_filename, 'r')
            self._blast_results = list(SeqIO.read(sequences_file, 'gb', IUPAC.ambiguous_dna))
            sequences_file.close()
        except Exception:
            print '\nFailed to load fetched sequences.'
    #end def
#end class