#!/usr/bin/python
# coding=utf-8

"""
Created on Mar 19, 2016

@author: Allis Tauri <allista@gmail.com>
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from BioUtils.SeqUtils import SeqLoader, get_indexes_of_genes
from BioUtils.NCBI import BlastCLI, BlastFilter


class Cluster(object):
    def __init__(self, genome, genes):
        self.genome = SeqLoader.load_file(genome)
        if not self.genome: raise RuntimeError('No genome loaded')
        self.indexes = [(i, get_indexes_of_genes(rec, genes)) for i, rec in enumerate(self.genome)]
        if len(self.indexes) > 1:
            raise NotImplementedError('Cluster spans across segments of the genome assembly.')
        self.contig  = self.genome[self.indexes[0][0]]
        self.indexes = self.indexes[0][1]
        #sanity check
        self.direction = np.sign(self.indexes[-1]-self.indexes[0])
        if not all((self.indexes[i+1]-self.indexes[i])*self.direction > 0 for i in xrange(len(self.indexes)-1)):
            raise NotImplementedError('Provided _genes do not follow each other in the genome')
        self.genes = [self.contig.features[i] for i in self.indexes]
        if len(self.genes) < 2: raise ValueError('Need at least 2 _genes to define a cluster')
        self.gaps = [self.genes[g+1].location.start-self.genes[g].location.end
                     if self.direction > 0 else
                     self.genes[g].location.start-self.genes[g+1].location.end
                     for g in xrange(len(genes)-1)]

    def __str__(self):
        genes = ['%s\ngap: %d' % (str(self.genes[g]), self.gaps[g]) for g in xrange(len(self.genes)-1)]
        genes.append(str(self.genes[-1]))
        return '\n'.join(genes)

    def __len__(self): return len(self.genes)

    def blast_genes(self, db, evalue=10, **kwargs):
        gresults = [None]*len(self)
        for g in self.genes:
            s = g.extract(self.contig)
            results = BlastCLI.blast_seq(s, db, evalue, **kwargs)
            print
            print 'Results: %s' % str(results)
            if not results:
                print 'No results'
                continue
            print g
            for rec in results:
                print 'Record: %s, %d alignments' % (str(rec), len(rec.alignments))
                for alignment in rec.alignments:
                    print 'Alignment: %s' % alignment
                    for hsp in alignment.hsps:
                        print 'HSP: %s' % hsp


class FuzzyValue(object):
    pvalue = 0.05

    def __init__(self, samples):
        self.samples = samples
        self.mean = np.mean(samples)
        self.SEM  = stats.sem(samples)
        self.min  = np.min(samples)
        self.max  = np.max(samples)

    def __eq__(self, other):
        t = stats.ttest_rel(self.samples, other.samples)
        return t[1] > self.pvalue

    def fits(self, value):
        return self.min*0.9 < value < self.max*1.1

    def __str__(self):
        return ('%f +/- %f (SEM)\nhist: %s'
                % (self.mean, self.SEM, stats.histogram(self.samples).count))


class ClusterModel(object):
    def __init__(self, clusters):
        self.ngens = len(clusters[0])
        if not all(len(c) == self.ngens for c in clusters[1:]):
            raise ValueError('Provided clusters have different number of _genes')
        self.lengths = [FuzzyValue([len(c.genes[i]) for c in clusters]) for i in xrange(self.ngens)]
        self.gaps = [FuzzyValue([c.gaps[i] for c in clusters]) for i in xrange(self.ngens-1)]

    def __str__(self):
        genes = ['gene: %s\ngap: %s' % (self.lengths[g], self.gaps[g]) for g in xrange(self.ngens-1)]
        genes.append('gene: %s' % str(self.lengths[-1]))
        return '\n'.join(genes)

    def cluster_fits(self, cluster):
        if len(cluster) != self.ngens: return False
        return (all(l.fits(len(g)) for l, g in zip(self.lengths, cluster.genes)) and
                all(g.fits(cg) for g, cg in zip(self.gaps, cluster.gaps)))


import os

if __name__ == '__main__':
    genome_dir = '/home/allis/Dropbox/Science/Микра/Thermococcus/sequence/GenBank/Thermococcus/'
    genes = {'Thermococcus_onnurineus_NA1-complete-genome.gb'
             : ['TON_1571', 'TON_1573'],
             'Thermococcus_barophilus_Ch5-complete.gb'
             : ['TBCH5v1_1758', 'TBCH5v1_1756'],
             'Thermococcus_barophilus_DT4-complete-genome.gb'
             : ['TBDT4v1_1587', 'TBDT4v1_1585'],
             'Thermococcus_sp._ES1.gb'
             : ['TES1_1488', 'TES1_1486'],
             'Thermococcus_gammatolerans_EJ3.gb'
             : ['TGAM_0057', 'TGAM_0055'],
             }

    clusters = [Cluster(os.path.join(genome_dir, filename), genes[filename]) for filename in genes]
    print '\n'.join(str(c) for c in clusters)
    model = ClusterModel(clusters)
    print model
    print
    print [model.cluster_fits(c) for c in clusters]

    clusters[0].blast_genes('nt', 10, remote=True)
