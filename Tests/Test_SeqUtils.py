'''
Created on Jan 15, 2016

@author: allis
'''
from BioUtils.Tools.Multiprocessing import MPMain
from BioUtils.Tools.Output import simple_timeit

class Test(MPMain):
    large_seqdb = u'/home/allis/Downloads/Local/Bioinformatics/SILVA/SILVA_123_SSURef_Nr99_tax_silva.fasta'
    #silva.gold.ng.fasta'#nogap.archaea.fasta' #SILVA_123_SSURef_Nr99_tax_silva.fasta'
    def _main(self):
        from BioUtils.SeqUtils import SeqView
        from BioUtils.Tools.Multiprocessing import parallelize_work
        
        with simple_timeit('load'):
            sv = SeqView()
            sv.load([self.large_seqdb])
            
        ssv = sv.subview(sv.keys()[:5])
        print ssv.keys()
        print ssv[3]
        print
        
        import cPickle as pickle
        ssv1 = pickle.loads(pickle.dumps(ssv, protocol=-1))
        print ssv1.keys()
        print ssv1[3]
        print
            
        def worker(r): return r.id    
        for numrecs in xrange(1000, len(sv), 1000):
            svs = sv[0:numrecs]
            with simple_timeit('sequential %d' % numrecs):
                res1 = [svs[i].id for i in range(numrecs)]
            with simple_timeit('parallel %d' % numrecs):
                res2 = parallelize_work(self.abort_event, 1, 1, worker, svs, copy_data=True)
            assert res1 == res2
            print '-'*80
        print 'Done'

if __name__ == '__main__':
    Test()
