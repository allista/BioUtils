# coding=utf-8

'''
Created on Dec 25, 2015

@author: Allis Tauri <allista@gmail.com>
'''

import sys, os

from Bio.Align import MultipleSeqAlignment
from Bio.SearchIO.HmmerIO import Hmmer3TextParser
from Bio.SeqFeature import FeatureLocation

from BioUtils.Tools.Multiprocessing import MultiprocessingBase, cpu_count
from BioUtils.Tools.Text import FilenameParser
from BioUtils.Tools.Misc import mktmp_name, run_cline
from BioUtils.Tools.Output import user_message
from BioUtils.SeqUtils import mktmp_fasta, Translator, get_indexes_of_all_genes
from BioUtils.HMMER3.Applications import HMMSearchCommandline, HMMBuildCommandline
from BioUtils.AlignmentUtils import AlignmentUtils
from BioUtils.Annotation import HSPAnnotator


class Hmmer(MultiprocessingBase, HSPAnnotator):
    
    def __init__(self, abort_event):
        super(Hmmer, self).__init__(abort_event)
        HSPAnnotator.__init__(self)
        
    @staticmethod
    def hmmbuild(alignment, outfile, name=None, **kwargs):
        unlink_file = False
        if isinstance(alignment, str): msafile = alignment
        elif isinstance(alignment, MultipleSeqAlignment):
            msafile = AlignmentUtils.mktmp(alignment)
            if not msafile: return False
            unlink_file = True
        else:
            print 'Alignment must be either a filename or an instance of MultipleSeqAlignment'
            return False
        if not name: name = FilenameParser.strip_ext(os.path.basename(outfile))
        ret = run_cline(HMMBuildCommandline(input=msafile, out=outfile,
                                            n=name, cpu=cpu_count, seed=0, **kwargs),
                        _msg = 'Unable to build HMM profile')
        if unlink_file: os.unlink(msafile)
        return ret

    @staticmethod
    def hmmsearch_recs(hmm, recs, **kwargs):
        recfile = mktmp_fasta(recs)
        hmm_out = mktmp_name('.hmm.txt')
        try:
            cline = HMMSearchCommandline(hmmfile=hmm, seqdb=recfile,
                                         o=hmm_out, cpu=cpu_count, seed=0,
                                         **kwargs)
            stdout, stderr = cline()
            print stdout
            if stderr: 
                sys.stderr.write(stderr)
                sys.stderr.flush()
            #parse hmmsearch results
            with open(hmm_out) as inp:
                return list(Hmmer3TextParser(inp))
        except Exception, e:
            print 'Error while running hmmsearch.'
            print e
            return None
        finally:
            os.unlink(recfile)
            os.unlink(hmm_out)

    annotation_type = 'hmmer'

    def hsp_score(self, hsp):
        return float(hsp.bitscore)

    def hsp2feature(self, name, group, location, hsp):
        feature = super(Hmmer, self).hsp2feature(name, group, location, hsp)
        feature.qualifiers['hmm_model'] = name
        feature.qualifiers['evalue'] = hsp.evalue
        feature.qualifiers['evalue_cond'] = hsp.evalue_cond
        feature.qualifiers['acc_average'] = hsp.acc_avg
        feature.qualifiers['bias'] = hsp.bias
        return feature

    def hmmsearch_genes(self, hmms, genome, table='Standard', decorate=False, **kwargs):
        #get _genes
        genes = get_indexes_of_all_genes(genome)
        if not genes: return None
        for gene_id, gi in enumerate(genes):
            genome.features[gi].qualifiers['feature_id'] = gi
            genome.features[gi].qualifiers['gene_id'] = gene_id
        #translate _genes
        with user_message('Translating _genes/CDS of %s' % genome.description, '\n'):
            translator = Translator(self._abort_event)
            translation = translator.translate_features(genome, genes, table)
        if not translation: return None
        if isinstance(hmms, str): hmms = [hmms]
        results = dict()
        for hmm in hmms:
            with user_message('Performing hmm search.'):
                hmm_results = self.hmmsearch_recs(hmm, translation, **kwargs)
            if not hmm_results: return None
            with user_message('Parsing search results...'):
                #get hit_ids of hmm matches
                hits = dict()
                for result in hmm_results:
                    for hit in result.iterhits():
                        hits[hit.id] = hit
                #get indexes of features where hmm hit
                hit_features = dict()
                for t in translation:
                    if t.id in hits:
                        fid = t.features[0].qualifiers.get('feature_id')
                        if fid is None: continue
                        hit_features[fid] = hits[t.id], t
                if hit_features: results.update(hit_features)
            #decorate genome
            if decorate:
                with user_message('Adding results as annotations...'):
                    hmm_name = os.path.basename(hmm)
                    for f in hit_features:
                        feature = genome.features[f]
                        for hsp in hit_features[f][0]:
                            if feature.strand == 1:
                                hmm_location = FeatureLocation(feature.location.start+hsp.hit_start*3,
                                                               feature.location.start+hsp.hit_end*3,
                                                               feature.strand)
                            else:
                                hmm_location = FeatureLocation(feature.location.end-hsp.hit_end*3,
                                                               feature.location.end-hsp.hit_start*3,
                                                               feature.strand)
                            hmm_feature = self.hsp2feature(hmm_name, 'HMM_annotations', hmm_location, hsp)
                            genome.features.append(hmm_feature)
        return results if results else None

    def hmmsearch_genome(self, hmms, genome, table='Standard', decorate=False, **kwargs):
        #translate _genes
        with user_message('Translating whole genome in 6 reading frames', '\n'):
            translator = Translator(self._abort_event)
            translation = translator.translate_six_frames(genome, table)
        if not translation: return None
        if isinstance(hmms, str): hmms = [hmms]
        results = []
        for hmm in hmms:
            with user_message('Performing hmm search.'):
                hmm_results = self.hmmsearch_recs(hmm, translation, **kwargs)
            if not any(len(r) for r in hmm_results): continue
            results += hmm_results
            #decorate genome
            if decorate:
                translation = dict((t.id, t) for t in translation)
                with user_message('Adding results as annotations...'):
                    hmm_name = os.path.basename(hmm)
                    glen = len(genome)
                    for frame in hmm_results:
                        for hit in frame:
                            frec = translation[hit.id]
                            start = frec.annotations['start']
                            strand = frec.annotations['strand']
                            for hsp in hit:
                                if strand == 1:
                                    hmm_location = FeatureLocation(start+hsp.hit_start*3,
                                                                   start+hsp.hit_end*3,
                                                                   strand)
                                else:
                                    hmm_location = FeatureLocation(glen-start-hsp.hit_end*3,
                                                                   glen-start-hsp.hit_start*3,
                                                                   strand)
                                hmm_feature = self.hsp2feature(hmm_name, 'HMM_annotations', hmm_location, hsp)
                                genome.features.append(hmm_feature)
        return results if results else None
