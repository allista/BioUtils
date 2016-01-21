#!/usr/bin/python
# coding=utf-8
#
# Copyright (C) 2015 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# degen_primer is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
Created on Dec 11, 2015

@author: Allis Tauri <allista@gmail.com>
'''

import os, re
from BioUtils.Tools.Output import isatty
from BioUtils.SeqUtils import load_files
from BioUtils.NCBI import BlastCLI

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC, NucleotideAlphabet, ProteinAlphabet

from BioUtils.Tools.Multiprocessing import MultiprocessingBase

class ClusterProject(MultiprocessingBase):
    
    def __init__(self, project, dirname, abort_event):
        super(ClusterProject, self).__init__(abort_event)
        self.proj = project
        self.dirname = dirname
        self.savefile = os.path.join(self.dirname, project+'.gb')
        self.genomes_files = []
        self.clusters = {}
        self.diagram = None
        self.fsets = {}
        self.genes = {}
        self.colors = {}
        self.crosslinks = {}
        try: os.mkdir(self.dirname)
        except: pass
        
    _illegal = re.compile(r'[ \-=/\(\)#]')
    _fidq = 'fid'
    
    @staticmethod
    def _feature_name(f, quals=('gene', 'locus_tag'), default='unknown'):
        for q in quals:
            if q in f.qualifiers:
                return ' '.join(f.qualifiers[q])
        return default

    def _process_features(self, rec):
        i = 0
        for f in rec.features:
            new_quals = dict()
            for q in f.qualifiers:
                qual = f.qualifiers[q]
                new_q = self._illegal.sub('_', q)
                new_quals[new_q] = qual
            f.qualifiers = new_quals
            if f.type == 'CDS':
                f.qualifiers[self._fidq] = [str(i)]
                i += 1
    
    def _extract_clusters(self, tag, qual='ugene_name'):
        tagre = re.compile(tag)
        clusters = {}
        for record in load_files(self._abort_event, self.genomes_files, 'gb'):
            for f in record.features:
                if qual in f.qualifiers:
                    q = ' '.join(f.qualifiers[qual])
                    if not tagre.match(q): continue
                    c = f.extract(record)
                    c.id = c.name = q
                    c.description = record.description
                    if c.seq.alphabet is not NucleotideAlphabet \
                    or c.seq.alphabet is not ProteinAlphabet:
                        c.seq.alphabet = IUPAC.IUPACAmbiguousDNA()
                    self._process_features(c)
                    clusters[c.id] = c
        return clusters
    
    def _extract_and_save(self, tag, qual='ugene_name', redo = False):
        clusters = None
        if not os.path.isfile(self.savefile) or redo:
            clusters = self._extract_clusters(tag, qual)
            SeqIO.write(clusters.values(), self.savefile, 'gb')
            print 'Clusters have been saved to: %s' % self.savefile
            return clusters
        else:
            clusters = dict((r.id, r) for r in SeqIO.parse(self.savefile, 'gb'))
            if not clusters:
                print 'No clusters have been found in "%s"' % self.savefile
                if isatty:
                    inp = input('Delete this file? [y/n]: ')
                    if inp.lower() == 'y': 
                        os.unlink(self.savefile)
                        print 'You can rerun the script now.'
                        return None
                print 'Delete this file manually and try again.'
                return None
        return clusters
    
    def _get_feature(self, cid, fid):
        if cid not in self.clusters: return None
        c = self.clusters[cid]
        for f in c.features:
            if self._fidq in f.qualifiers \
            and f.qualifiers[self._fidq][0] == str(fid):
                return f
    
    @MultiprocessingBase.data_mapper_method
    def _blast_feature(self, f, c1, c2, features1, features2, evalue, max_rlen):
        results = BlastCLI.s2s_blast(f.extract(c1), c2, evalue, command='blastn', task='blastn')
        hsps = BlastCLI.all_hsps(results, max_rlen)
        if not hsps: return [(None, None, None)]
        f1 = []
        f2 = []
        col = []
        for hsp in hsps:
            col.append(colors.linearlyInterpolatedColor(colors.Color(1,1,1,0.2), colors.Color(0,0,0,0.2), 
                                                        0, 1, float(hsp.identities)/hsp.align_length))
            f1.append(SeqFeature(FeatureLocation(f.location.start+hsp.query_start, 
                                                 f.location.start+hsp.query_start+hsp.align_length, strand=0)))
            f2.append(SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_start+hsp.align_length, strand=0)))
        return zip(f1, f2, col)
    
    @MultiprocessingBase.results_assembler_methd
    def _compose_crosslink(self, index, result, features1, features2):
        for f1, f2, col in result:
            if f1 is None: continue
            tf1 = features1.add_feature(f1, color=col, border=col)
            tf2 = features2.add_feature(f2, color=col, border=col)
            self.diagram.cross_track_links.append(CrossLink(tf1, tf2, col, col))
    
    def _compose_crosslinks(self, cluster_ids, evalue, max_rlen):
        if not self.diagram or not self.fsets: return
        print 'Finding crosslinks between genes in the cluster.'
        num_cids = len(cluster_ids)
        for ci, cid1 in enumerate(cluster_ids):
            if ci >= num_cids-1: break
            features1 = self.fsets[cid1]
            c1 = self.clusters[cid1]
            cid2 = cluster_ids[ci+1]
            features2 = self.fsets[cid2]
            c2 = self.clusters[cid2]
            print '%s vs %s' % (cid1, cid2)
            work = self.Work()
            work.prepare_jobs(self._blast_feature, 
                              [f for f in c1.features if f.type == 'CDS'], None, 
                              c1, c2, features1, features2, evalue, max_rlen)
            work.set_assembler(self._compose_crosslink, features1, features2)
            self.start_work(work)
            self.wait(work)
            
    def get_clusters(self, genomes_dir, tag, force = False, print_ids = False):
        print 'Extracting clusters from provided GenBank files.'
        self.genomes_files = [os.path.join(genomes_dir, f) 
                              for f in os.listdir(genomes_dir) if f.endswith('.gb')]
        self.clusters = self._extract_and_save(tag, redo=force)
        if not self.clusters: return False
        if print_ids: print self.clusters.keys()
        return True
        
    def draw_clusters(self, cluster_ids=None, evalue=0.0001, max_rlen=0.01, border=True, pagesize='A4', add_crosslinks=True):
        print 'Creating cluster diagram.'
        #create diagram
        self.diagram = GenomeDiagram.Diagram(self.proj)
        #add tracks
        max_len = 0
        if cluster_ids is None:
            cluster_ids = sorted(self.clusters.keys())
        for cid in cluster_ids:
            cluster = self.clusters[cid]
            clen = len(cluster)
            if clen > max_len: max_len = clen
            track = self.diagram.new_track(1,
                                      name=cluster.description,
                                      greytrack=True, height=0.4,
                                      greytrack_fontcolor=colors.black,
                                      greytrack_labels=1,
                                      scale=False,
                                      start=0, end=clen)
            self.fsets[cid] = track.new_set()
        #add crosslink features
        if add_crosslinks:
            self._compose_crosslinks(cluster_ids, evalue, max_rlen)
        #add CDS-es
        for cid in cluster_ids:
            cluster = self.clusters[cid]
            genes = self.genes.get(cluster.id, [])
            cols = self.colors.get(cluster.id, [])
            i = 1
            for f in cluster.features:
                if f.type != 'CDS': continue
                fname = genes[i-1] if i <= len(genes) else 'CDS_%d' % i
                fcol = cols[i-1] if i <= len(cols) else colors.grey
                self.fsets[cid].add_feature(f, sigil="BIGARROW",
                                            color=fcol, 
                                            border = colors.black if border else fcol,
                                            name = fname, label=True,
                                            label_position="middle",
                                            label_size = 8, label_angle=45)
                i += 1
        self.diagram.draw(format="linear", pagesize=pagesize, fragments=1,
                     start=0, end=max_len)
        for ptype in ('PDF', 'EPS', 'SVG'):
            self.diagram.write(os.path.join(self.dirname, '%s.%s' %
                                            (self.proj, ptype.lower())), ptype)
#end def

def grey(val): return colors.Color(val, val, val)

def FC_cluster(abort_event, color = True, force = False, print_ids = False, add_crosslinks=True):
    proj = ClusterProject('Formate Clusters', u'/home/allis/Dropbox/Science/Микра/Thermococcus/Figures/FC', abort_event)
    if not proj.get_clusters(genomes_dir, '.+\-FC\-full', force=force, print_ids=print_ids): return
    
    #formate cluster specifics
    other_genes = ('fdhA', '4Fe-4S', 
                  'mbhH', 'mbhH\'', 'mbhH\'\'',
                  'mbhM', 'mbh K+L', 'mbhN', 
                  'mbhJ', 'mbhX', 'focA', 
                  'mbhB', 'mbhC', 'mbhD', 'mbh E+F', 'mbhG', 'mbhA', 'mbhH\'\'\'')
    
    Ch5_genes = list(other_genes)
    Ch5_genes.insert(4, 'mbhH\'')
    
    if color:
        other_colors = [colors.Color(1,69.0/255)]*2 + \
                     [colors.red]*6 + \
                     [colors.Color(1,185.0/255), colors.darkgreen, colors.darkmagenta] + \
                     [colors.aqua]*7
                     
        Ch5_colors = list(other_colors)
        Ch5_colors.insert(4, colors.red)
                     
        
    else:
        other_colors = [grey(0.2)]*2 + \
                     [grey(0.5)]*6 + \
                     [grey(0.3), grey(0.6), grey(0.25)] + \
                     [grey(0.95)]*7
                     
        Ch5_colors = list(other_colors)
        Ch5_colors.insert(4, grey(0.5))
        
                 
    cluster_ids = ('TBCH5-FC-full',
                   'TONA1-FC-full',
                   'TGAM-FC-full',
                   'TES1-FC-full',
                   'TBDT4-FC-full',
                   'THDS1-FC-full',
                   )
    
    for cid in cluster_ids:
        if cid == 'TBCH5-FC-full':
            proj.genes[cid] = Ch5_genes
            proj.colors[cid] = Ch5_colors
        else:
            proj.genes[cid] = other_genes
            proj.colors[cid] = other_colors
    
    proj.draw_clusters(cluster_ids, add_crosslinks=add_crosslinks)
#end

def CO_cluster(abort_event, color = True, force = False, print_ids = False, add_crosslinks=True):
    proj = ClusterProject('CO Clusters', u'/home/allis/Dropbox/Science/Микра/Thermococcus/Figures/CO', abort_event)
    if not proj.get_clusters(genomes_dir, '.+\-COC\-full', force=force, print_ids=print_ids): return
    
    #CO cluster specifics
    other_genes = ('corQ', 'corR', 'cooF', 'cooS', 'cooC',
                   'unknown',
                   'mbhH', 'mbhM', 'mbh K+L', 'mbhN', 'mbhJ',
                   'mbhB', 'mbhC', 'mbhD', 'mbh E+F', 'mbhG', 'mbhA', 'mbhH\'\'\'')
    
    AM4_genes = list(other_genes)
    AM4_genes.insert(5, 'unknown')
    
    if color: #TODO: change colors for CO cluster
        other_colors = [colors.Color(1,69.0/255)]*2 + \
                     [colors.red]*6 + \
                     [colors.Color(1,185.0/255), colors.darkgreen, colors.darkmagenta] + \
                     [colors.aqua]*7
                     
        AM4_colors = list(other_colors)
        AM4_colors.insert(2, colors.red)
                     
        
    else:
        other_colors = [grey(0.2)]*2 + \
                     [grey(0.3)]*3 + \
                     [grey(0.1)]+ \
                     [grey(0.5)]*5 + \
                     [grey(0.95)]*7
                     
        AM4_colors = list(other_colors)
        AM4_colors.insert(5, grey(0.1))
        
                 
    cluster_ids = ('TBCH5-COC-full',
                   'TBMP-COC-full',
                   'TAM4-COC-full',
                   'TONA1-COC-full',
                   )
    
    for cid in cluster_ids:
        if cid == 'TAM4-COC-full':
            proj.genes[cid] = AM4_genes
            proj.colors[cid] = AM4_colors
        else:
            proj.genes[cid] = other_genes
            proj.colors[cid] = other_colors
    
    proj.draw_clusters(cluster_ids, add_crosslinks=add_crosslinks)
#end

if __name__ == '__main__':
    from multiprocessing import Event
    abort_event = Event()
    genomes_dir = u'/home/allis/Dropbox/Science/Микра/Thermococcus/sequence/GenBank/Thermococcus'
    FC_cluster(abort_event, color=False, force=False, add_crosslinks=True)
#    CO_cluster(color=False, force=False, add_crosslinks=False)
    print 'Done'