# coding=utf-8

'''
Created on Mar 17, 2016

@author: Allis Tauri <allista@gmail.com>
'''

from time import time, sleep
from datetime import timedelta

from Bio import SeqIO
from Bio import Entrez

from BioUtils.Tools.Misc import retry

class BatchEntrez(object):
    #defaults
    RETRIES    = 3
    PAUSE_EACH = 100
    BATCH      = 20
    PAUSE      = 60
    
    def __init__(self, email):
        self.email       = email
        self._start_time = -1
    #end def
    
    def _get_records(self, query, database, **efetch_kwargs):
        results = retry(lambda : Entrez.read(Entrez.esearch(db=database, term=query, usehistory="y")),
                        'Unable to get Entrez ID for query: %s' % query, self.RETRIES)
        if not results['IdList']:
            print('NCBI returned no result for query: %s' % query) 
            return []
        webenv    = results['WebEnv']
        query_key = results['QueryKey']
        #fetch genbank data for the received IDs 
        num_results = len(results['IdList'])
        print('Downloading data...')
        data = retry(lambda : Entrez.efetch(db=database, rettype="gb", retmode="text",
                                            retstart=0, retmax=num_results,
                                            webenv=webenv, query_key=query_key, **efetch_kwargs),
                     'Unable to download data for IDs: %s' % str(results['IdList']), self.RETRIES)
        #parse received data
        try: records = list(SeqIO.parse(data, 'gb'))
        except Exception as e:
            print 'Unable to parse fetched data as SeqRecords'
            print e
            return []
        finally: data.close()
        print('Done. Elapsed time: %s\n' % timedelta(seconds=time()-self._start_time))
        return records
    
    def _get_records_for_terms(self, terms, start, stop, database, **efetch_kwargs):
        query = ' or '.join(terms[i] for i in xrange(start, stop))
        num_terms = len(terms)
        print('[%3.1f%%] performing query for terms %d-%d of %d' 
              % (float(stop)/num_terms*100, start+1, stop, num_terms))
        return self._get_records(query, database, **efetch_kwargs)
    #end def
    
    def get_records(self, query, database, **efetch_kwargs):
        self._start_time = time()
        Entrez.email = self.email
        Entrez.tool = 'BioUtils.BatchEntrez.get_records'
        return self._get_records(query, database, **efetch_kwargs)
    
    def get_records_for_terms(self, terms, database, **efetch_kwargs):
        self._start_time = time()
        #check number of queries and warn the user
        num_terms   = len(terms)
        num_queries = num_terms/self.BATCH
        num_pauses  = 0; pause_time = 0 
        if num_queries > self.PAUSE_EACH:
            num_pauses = num_queries/self.PAUSE_EACH
            pause_time = num_pauses * self.PAUSE
            self.PAUSE_EACH = num_queries/(num_pauses+1)+1
            print('WARNING: %d separate Entrez queries will be made.\n'
                  'To comply with NCBI rules the queries will be made\n'
                  'in series of %d with %d sec pause in between.\n' 
                  % (num_queries, self.PAUSE_EACH, self.PAUSE))
            print('Total pause time will be:\n%s\n' % timedelta(seconds=pause_time))
        query_time = num_queries * 4/3.0
        if query_time > 5:
            print('No more than 3 requests per second is allowed by NCBI,\n'
                  'so *minimum* time spend for your query will be:\n%s\n' % timedelta(seconds=query_time))
        if pause_time > 0 and query_time > 5:
            print('Total *minimum* estimated time:\n%s\n' % timedelta(seconds=pause_time+query_time))
            print('Note, that depending on the load of NCBI servers it\n'
                  'may take several times as much.\n')
        #setup Entrez engine
        Entrez.email = self.email
        Entrez.tool = 'BioUtils.BatchEntrez.get_records_for_terms'
        #perform queries in batches
        pause_num = self.PAUSE_EACH
        records = []
        for i in xrange(0, num_terms, self.BATCH):
            if i/self.BATCH > pause_num:
                print('Pausing for %d seconds...\n' % self.PAUSE)
                sleep(self.PAUSE)
                pause_num += self.PAUSE_EACH
            try: records.extend(self._get_records_for_terms(terms, i, min(i+self.BATCH, num_terms), database, **efetch_kwargs))
            except RuntimeError as e: 
                print(e)
                continue
        print('Done.\nTotal elapsed time: %s\n' % timedelta(seconds=time()-self._start_time))
        return records
    #end def
#end class
