import os
import re
import shlex
from bisect import bisect_left
from collections import Counter
from glob import glob
from itertools import chain
from reportlab.lib import colors
from ConfigParser import SafeConfigParser, Error as ConfigError

from Bio.Alphabet import NucleotideAlphabet
from Bio.Alphabet import generic_dna
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation

from BioUtils.Annotation import AnnotatorBase, HSPAnnotator
from BioUtils.NCBI import BlastCLI
from BioUtils.SeqUtils import SeqView, filename_for_record, all_CDS, rec_name_with_id, backup_write
from BioUtils.SeqUtils import Translator, copy_attrs
from BioUtils.Tools.Multiprocessing import MultiprocessingBase, ordered_results_assembler
from BioUtils.Tools.Multiprocessing import ordered_shelved_results_assembler
from BioUtils.Tools.Output import ProgressCounter
from BioUtils.Tools.tmpStorage import shelf_result


class TaggedGene(object):
    class Tag(object):
        def __init__(self, name, coverage):
            self.name = name
            self.coverage = coverage

        def better(self, other):
            return self.coverage > other.coverage

        def update(self, coverage):
            self.coverage = max(self.coverage, coverage)

        def __str__(self):
            return '%s %.1f' % (self.name, self.coverage)

        def __repr__(self): return str(self)

    def __init__(self, index, gene, tag_feature, intersection):
        self.gi = index
        self.gene = gene
        self._glen = float(len(gene))
        self.tag = None
        self.update_tag(tag_feature, intersection)

    @property
    def start(self):
        return self.gene.location.start

    @property
    def end(self):
        return self.gene.location.end

    @property
    def strand(self):
        return self.gene.location.strand

    def update_tag(self, tag_feature, intersection):
        name = AnnotatorBase.get_tag(tag_feature)
        if not name:
            raise ValueError('tag_feature should contain %s qualifier' % AnnotatorBase.tag_qualifier)
        coverage = HSPAnnotator.get_score(tag_feature) * intersection/self._glen
        tag = self.Tag(name, coverage)
        if self.tag is None or tag.better(self.tag):
            self.tag = tag

    def update(self, other):
        if self.gi != other.gi: return
        if self.tag is None or other.tag.better(self.tag):
            self.tag = other.tag


    def __str__(self): return str(self.tag)

    def __repr__(self): return str(self)


class GenesTagger(object):
    def __init__(self, gene_features, parent_name):
        self.pname = parent_name
        #check for compound locations
        new_genes = []
        remove = set()
        gene_features = set(gene_features)
        for f in gene_features:
            if isinstance(f.location, CompoundLocation):
                for part in f.location.parts:
                    new_f = SeqFeature(part,f.type,id=f.id,qualifiers=f.qualifiers)
                    new_genes.append(new_f)
                remove.add(f)
        if new_genes:
            gene_features -= remove
            gene_features.update(new_genes)
        self._genes = sorted(gene_features, key=lambda a: a.location.end)
        self._ends = [f.location.end for f in self._genes]
        self._len = len(self._genes)
        self._last = self._len-1

    def __getitem__(self, item):
        return self._genes[item]

    def __len__(self): return self._len

    @staticmethod
    def _intersection(a, b):
        s = max(a.location.start, b.location.start)
        e = min(a.location.end, b.location.end)
        return e-s if s < e else 0

    def _insert_new_CDS(self, i, f):
        # print 'WARNING: %s probably have missing CDS annotation at: %s' % (self.pname, str(f.location))
        f.type = 'CDS'
        self._genes.insert(i, f)
        self._ends.insert(i, f.location.end)
        self._len += 1
        self._last += 1

    def tagged_genes(self, tag, min_coverage=0, same_strand=False):
        i = bisect_left(self._ends, tag.location.start)
        if i == self._last:
            self._insert_new_CDS(i, tag)
            return [TaggedGene(i, tag, tag, len(tag))]
        intersection = []
        for gi, g in enumerate(self._genes[i:]):
            if same_strand and g.location.strand != tag.location.strand: continue
            _i = self._intersection(tag, g)
            if _i: intersection.append((gi+i, g, _i))
        if not intersection:
            tag.type = 'CDS'
            self._insert_new_CDS(i, tag)
            return [TaggedGene(i, tag, tag, len(tag))]
        # print '\n'.join('%d %s\nintersection: %f > %f' % (gi, tag, float(_i)/len(g), min_coverage)
        #                 for gi, g, _i in intersection)#debug
        return [TaggedGene(gi, g, tag, _i)
                for gi, g, _i in intersection
                if float(_i)/len(g) > min_coverage]


class GeneCluster(object):
    def __init__(self, cluster_def):
        self.name = cluster_def.name
        self.max_gap = cluster_def.max_insertions
        self._genes = dict()
        self.max_gi = 0

    def __len__(self): return len(self._genes)

    def add(self, gene):
        if self._genes and gene.gi-self.max_gi-1 > self.max_gap: return False
        g = self._genes.get(gene.gi)
        if not g:
            self._genes[gene.gi] = gene
            self.max_gi = max(gene.gi, self.max_gi)
        else: g.update(gene)
        return True

    @property
    def gi_set(self):
        return set(self._genes.keys())

    @property
    def gis(self):
        return sorted(self._genes.keys())

    @property
    def genes(self):
        return [self._genes[gi] for gi in self.gis]

    @property
    def tags(self):
        return [self._genes[gi].tag.name for gi in self.gis]

    @property
    def location(self):
        if not self._genes: return None
        start = -1; end = -1
        strand = Counter()
        for g in self._genes.values():
            s = g.start; e = g.end
            if start < 0 or s < start: start = s
            if e > end: end = e
            strand[g.strand] += 1
        return FeatureLocation(start, end, strand.most_common(1)[0][0])

    def issubcluster(self, other):
        return self.gi_set.issubset(other.gi_set)

    def supercluster(self, other):
        sset = self.gi_set
        oset = other.gi_set
        if sset.issubset(oset): return other
        if oset.issubset(sset): return self
        return None

    def __str__(self):
        loc = self.location
        genes = self.genes
        if loc.strand == -1: genes = genes[::-1]
        return '%s: [%s]' % (str(self.location), ', '.join(str(g) for g in genes))

    def __repr__(self): return str(self)


class ClusterString(object):
    class _Tag(object):
        def __init__(self, tags):
            if not tags:
                raise  ValueError('tags should contain at least one tag')
            self.single = len(tags) == 1
            self.value = tags[0] if self.single else set(tags)

        def match(self, tag):
            if self.single:
                return tag == self.value
            return tag in self.value

        def __eq__(self, other):
            return self.match(other)

        def __str__(self):
            if self.single:
                return self.value
            return '|'.join(self.value)

    def __init__(self, string):
        self.tags = [self._Tag(t) for t in (t.strip().split('|') for t in string.split())]

    def missing(self, tags):
        self_set = set(self.tags)
        for other in tags:
            remove = None
            for tag in self_set:
                if tag.match(other):
                    remove = tag
                    break
            if remove:
                self_set.remove(remove)
                if not self_set: break
        return self_set

    def match(self, tags):
        return not self.missing(tags)

    def match_strict(self, tags):
        start = self.tags.index(tags[0])
        if start+len(tags) > len(self.tags): return False
        prev = start
        for t in tags[1:]:
            next = self.tags.index(t)
            if next < prev: return False
            prev = next
        return True

    def __contains__(self, item):
        return any(t.match(item) for t in self.tags)

    def __iter__(self):
        return self.tags.__iter__()

    def __str__(self):
        return ', '.join(str(t) for t in self.tags)


class ClusterDef(object):
    def __init__(self, name, full, core='', strict_order=False, max_insertions=0, min_coverage=0):
        self.name = name
        self.full_cluster = ClusterString(full)
        self.core_cluster = ClusterString(core) if core else self.full_cluster
        self.strict_order = strict_order
        self.max_insertions = max_insertions
        self.min_coverage = min_coverage

    def __contains__(self, item):
        return item in self.full_cluster

    def match(self, cluster):
        ctags = cluster.tags
        if self.core_cluster.missing(ctags): return False
        if self.strict_order:
            strand = cluster.location.strand
            if strand == -1: ctags = ctags[::-1]
            if not self.full_cluster.match_strict(ctags): return False
        return True


class ClusterFinder(MultiprocessingBase, AnnotatorBase):
    annotation_type = 'cluster'
    genes_qualifier = 'cluster_genes'

    def __init__(self, abort_event):
        super(ClusterFinder, self).__init__(abort_event)
        AnnotatorBase.__init__(self)
        self.clusters = []

    def load(self, cluster_files):
        parser = SafeConfigParser(defaults=dict(core='', strict_order=None, max_insertions='0', min_coverage='0.1'))
        try: parser.read(cluster_files)
        except ConfigError as e:
            print 'Error while parsing cluster configuration:\n%s' % str(e)
            return False
        self.clusters = []
        for name in parser.sections():
            try: full = parser.get(name, 'full')
            except ConfigError:
                print 'No "full" value in %s cluster definition.' % name
                continue
            try:
                core = parser.get(name, 'core')
                strict_order = parser.get(name, 'strict_order') is not None
                max_insertions = parser.getint(name, 'max_insertions')
                min_coverage = parser.getfloat(name, 'min_coverage')
            except (ConfigError, ValueError) as e:
                print '%s cluster definition is malformed:\n%s' % (name, str(e))
                continue
            self.clusters.append(ClusterDef(name, full, core, strict_order, max_insertions, min_coverage))
        if not self.clusters:
            print 'No cluster definitions were loaded from provided files.'
            return False
        return True

    @MultiprocessingBase.data_mapper_method
    def _process_genome(self, gi, db, output_dir):
        genome = db[gi]
        # sort features
        tags = sorted((f for f in genome.features if AnnotatorBase.tag_qualifier in f.qualifiers), key=lambda a: a.location.end)
        tagger = GenesTagger(all_CDS(genome), rec_name_with_id(genome))
        # search for clusters
        found = []
        report = []
        for cluster in self.clusters:
            candidate = GeneCluster(cluster)
            candidates = [candidate]
            for f in tags:
                tag = f.qualifiers.get(BlastCLI.tag_qualifier)
                if not tag: tag = f.qualifiers.get('gene')
                if not tag: continue
                tag = tag[0]
                if tag in cluster:
                    for tg in tagger.tagged_genes(f, cluster.min_coverage):
                        if candidate.add(tg): continue
                        candidate = GeneCluster(cluster)
                        candidate.add(tg)
                        candidates.append(candidate)
            # check cluster candidates
            confirmed = []
            for c in candidates:
                if not cluster.match(c): continue
                is_subcluster = False
                for fi in range(len(found)-1,-1,-1):
                    f = found[fi]
                    if c.issubcluster(f):
                        is_subcluster = True
                        break
                    if f.issubcluster(c):
                        del found[fi]
                if is_subcluster: continue
                confirmed.append(c)
            if confirmed: found += confirmed
        # save if found anything
        if found:
            for c in found:
                location = c.location
                cgenes = c.genes
                if location.strand == -1: cgenes = cgenes[::-1]
                gene_names = []
                last_gi = -1
                for g in cgenes:
                    gap = abs(g.gi-last_gi)
                    if last_gi >= 0 and gap > 1:
                        gene_names.extend(['NONE']*(gap-1))
                    gene_names.append(g.tag.name)
                    last_gi = g.gi
                feature = self.annotate_location(c.name, 'putative_clusters', c.location)
                feature.qualifiers[self.genes_qualifier] = ' '.join(gene_names)
                genome.features.append(feature)
                report.append('%s %s' % (str(location), c.name))
            fname = filename_for_record(genome, 'clusters.gb')
            fname = os.path.join(output_dir, fname)
            genome.seq.alphabet = generic_dna
            backup_write(genome, fname)
        return report, rec_name_with_id(genome)

    def find_clusters(self, files, output_dir=''):
        if not self.clusters: return False
        # load genomes
        genomes = SeqView.safe_load(files)
        if not genomes: return False
        # create dest dir if needed
        if output_dir and not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        # process files
        results = {}
        glen = len(genomes)
        with ProgressCounter('Searching for %d clusters in %d sequence%s:' %
                                     (len(self.clusters), glen, 's' if glen>1 else ''), glen) as prg:
            work = self.Work()
            work.start_work(self._process_genome, genomes.keys(), None, genomes, output_dir)
            work.assemble(ordered_results_assembler, results, prg)
            if not work.wait(): return False
        # check results
        found = False
        for gi in results:
            res, gname = results[gi]
            if not res: continue
            found = True
            print('%s:\n\t%s\n' % (gname, '\n\t'.join(res)))
        if not found:
            print 'No putative clusters were found.'
        return True


class _DiagramCluster(object):
    _illegal = re.compile(r'[ \-=/()#]')
    _trash_info = \
    [
        re.compile(r',?\s*whole genome shotgun sequence\.?'),
        re.compile(r',?\s*complete genome\.?'),
        re.compile(r'\s*genomic scaffold'),
        re.compile(r'\s*[Ss]caffold\w+\b'),
        re.compile(r'\s*strain\b'),
    ]

    @classmethod
    def _cleanup_sequence(cls, rec):
        #cleanup description
        desc = rec.description
        for exp in cls._trash_info:
            desc = exp.sub('', desc)
        rec.description = desc
        #cleanup features
        for f in rec.features:
            new_quals = dict()
            for q in f.qualifiers:
                qual = f.qualifiers[q]
                new_q = cls._illegal.sub('_', q)
                new_quals[new_q] = qual
            f.qualifiers = new_quals

    def __init__(self, seq, loc):
        self.sequence = seq
        self.genes = []
        self.marked = False
        self.location = loc
        if self.sequence.seq.alphabet is not NucleotideAlphabet:
            self.sequence.seq.alphabet = generic_dna
        self._cleanup_sequence(self.sequence)

    def __len__(self): return len(self.sequence)

    def set_genes(self, qualifier):
        if not qualifier: self.genes = []
        else: self.genes = [n.strip() for n in ' '.join(qualifier).split()]

    @property
    def label(self):
        strand = ''
        if self.location.strand == 1: strand = '+'
        elif self.location.strand == -1: strand = '-'
        return '%s, %s [%d..%d]%s' % (self.sequence.description, self.sequence.id,
                                       self.location.start+1, self.location.end, strand)

    @property
    def CDS(self):
        return all_CDS(self.sequence)


class ClusterDiagram(MultiprocessingBase):
    def __init__(self, abort_event):
        super(ClusterDiagram, self).__init__(abort_event)
        self.name = 'ClusterDiagram'
        self.project = None
        self.files = []
        #general options
        self.tags = []
        self.order = []
        self.marks = set()
        self.colors = {}
        # display options (also in general section)
        self.pagesize = 'A4'
        self.default_color = colors.grey
        self.no_border = False
        self.no_color = False
        #crosslinks options
        self.add_crosslinks = True
        self.min_identity = 0
        self.min_length = 0
        self.evalue = 10
        #label options
        self.name_size = 9
        self.gene_size = 10
        self.unknown_gene_size = 7
        self.gene_angle = 15
        #results
        self.clusters = []
        self.fsets = []
        self.crosslinks = []
        self.diagram = None

    def _get_option(self, section, option, default, conv=None):
        try:
            val = self.project.get(section, option)
            if conv is not None: val = conv(val)
        except (ConfigError, ValueError): val = default
        return val

    def load(self, project_file, overrides=None): #TODO: implement overrides
        self.project = SafeConfigParser()
        try: self.project.read(project_file)
        except ConfigError as e:
            print 'Unable to load %s:\n%s' % (project_file, str(e))
            return False
        # parse general section
        if not self.project.has_section('general'):
            print 'Malformed project file.\nNo "general" section in %s' % project_file
            return False
        self.tags = shlex.split(self._get_option('general', 'tags', ''))
        if not self.tags:
            print 'No "tags" were value was found in "general" section in %s' % project_file
            return False
        self.name = project_file
        self.order = shlex.split(self._get_option('general', 'order', ''))
        self.marks = set(shlex.split(self._get_option('general', 'marks', '')))
        self.pagesize = self._get_option('general', 'page_size', 'A4')
        self.no_border = self._get_option('general', 'no_border', None) is not None
        self.no_color = self._get_option('general', 'no_color', None) is not None
        self.default_color = self._get_option('general', 'default_color', '')
        if self.default_color:
            self.default_color = colors.HexColor(self.default_color)
        else: self.default_color = colors.grey
        # files section
        self.files = []
        if self.project.has_section('files'):
            dir = self._get_option('files', 'dir', '')
            files = shlex.split(self._get_option('files', 'files', ''))
            if files:
                self.files = list(chain.from_iterable(glob(os.path.join(dir, f)) for f in files))
        # crosslinks options
        self.add_crosslinks = self.project.has_section('crosslinks')
        if self.add_crosslinks:
            self.min_identity = self._get_option('crosslinks', 'min_identity', 0, float)
            self.min_length = self._get_option('crosslinks', 'min_length', 0, float)
            self.evalue = self._get_option('crosslinks', 'evalue', 10, float)
        # colors if have/needed
        if not self.no_color:
            if self.project.has_section('colors'):
                self.colors = dict((gene, colors.HexColor(col))
                                   for gene, col in self.project.items('colors'))
            #else, colors will be generated automatically after the clusters are extracted
        # label options
        self.name_size = max(6, self._get_option('labels', 'name_size', 10, int))
        self.gene_size = max(6, self._get_option('labels', 'gene_size', 10, int))
        self.unknown_gene_size = max(6, self.gene_size - 2)
        self.gene_angle = self._get_option('labels', 'gene_angle', 15, float)
        return True

    @staticmethod
    def _feature_name(f, quals=('gene', 'locus_tag'), default='unknown'):
        for q in quals:
            val = f.qualifiers.get(q)
            if not val: continue
            return ' '.join(val)
        return default


    @MultiprocessingBase.data_mapper_method
    @shelf_result
    def _process_genome(self, gi, recs):
        record = recs[gi]
        clusters = []
        for tag in self.tags:
            tag_clusters = []
            for f in record.features:
                # cluster annotation
                q = f.qualifiers.get(ClusterFinder.tag_qualifier)
                if not q: continue
                if tag != ' '.join(q): continue
                cluster = f.extract(record)
                copy_attrs(record, cluster)
                if tag_clusters: cluster.description += ' cluster %d' % (len(tag_clusters)+1)
                cluster = _DiagramCluster(cluster, f.location)
                cluster.set_genes(f.qualifiers.get(ClusterFinder.genes_qualifier, []))
                cluster.marked = record.id in self.marks
                tag_clusters.append(cluster)
            if tag_clusters: clusters += tag_clusters
        return clusters

    def _generate_gene_colors(self):
        if not self.clusters: return
        full_gene_set = set()
        for c in self.clusters:
            full_gene_set.update(c.genes)
        if 'NONE' in full_gene_set:
            full_gene_set.remove('NONE')
        ngenes = float(len(full_gene_set))-1
        middle = ngenes/2.0
        self.colors = {}
        for i, gene in enumerate(sorted(full_gene_set)):
            t = i/ngenes
            if i < middle:
                c = colors.linearlyInterpolatedColor(colors.Color(1, 0, 0, 1), colors.Color(0, 1, 0, 1), 0, 1, t*2)
            else:
                c = colors.linearlyInterpolatedColor(colors.Color(0, 0.9, 0.1, 1), colors.Color(0, 0, 1, 1), 0, 1, t*2-1)
            self.colors[gene] = c

    def extract_clusters(self):
        self.clusters = []
        genomes = SeqView.safe_load(self.files)
        if not genomes: return False
        glen = len(genomes)
        self.clusters = [None]*glen
        if self.order: self.order = [oid for oid in self.order if oid in genomes.keys()]
        with ProgressCounter('Extracting clusters from provided genomes:', glen) as prg:
            work = self.Work()
            work.start_work(self._process_genome, self.order or genomes.keys(), None, genomes)
            work.assemble(ordered_shelved_results_assembler, self.clusters, prg)
            if not work.wait(): return False
        self.clusters = list(chain.from_iterable(c for c in self.clusters if c))
        #generate gene colors if needed
        if not self.no_color and not self.colors:
            self._generate_gene_colors()
        return bool(self.clusters)

    @MultiprocessingBase.data_mapper_method
    def _blast_feature(self, f, c1, c2):
        trans = Translator(self._abort_event)
        cds = trans.translate(f.extract(c1), 11)
        sixframes = trans.translate_six_frames_single(c2, 11)
        if not sixframes: return [(None, None, None)]
        results = []
        for frame in sixframes:
            res = BlastCLI.s2s_blast(cds, frame, self.evalue, command='blastp', task='blastp')
            if res: results.extend(res)
        hsps = BlastCLI.all_hsps(results)
        if not hsps: return [(None, None, None)]
        f1 = []
        f2 = []
        col = []
        fname = self._feature_name(f, default='CDS')
        cds_len = len(cds)
        min_len = len(cds) * self.min_length
        for hsp in hsps:
            if hsp.align_length < min_len: continue
            if hsp.identities / float(hsp.align_length) < self.min_identity: continue
            color_t = (float(hsp.identities) / hsp.align_length)
            print '%s %s: %5.1f%% (%5.1f%%)' % (c1.description, fname, color_t * 100, float(hsp.identities) / cds_len * 100)
            col.append(colors.linearlyInterpolatedColor(colors.Color(0, 0, 1, 0.2), colors.Color(0, 1, 0, 0.2),
                                                        0.2, 1, color_t))
            qstart = (hsp.query_start - 1) * 3
            qend = qstart + hsp.align_length * 3
            sstart = (hsp.sbjct_start - 1) * 3
            send = sstart + hsp.align_length * 3
            f1.append(
                SeqFeature(FeatureLocation(f.location.start + qstart, f.location.start + qend, strand=hsp.strand[0])))
            f2.append(SeqFeature(FeatureLocation(sstart, send, strand=hsp.strand[1])))
        return zip(f1, f2, col)

    @MultiprocessingBase.results_assembler_method
    def _compose_crosslink(self, index, result, features1, features2):
        for f1, f2, col in result:
            if f1 is None: continue
            if self.no_color: col = colors.color2bw(col)
            tf1 = features1.add_feature(f1, color=col, border=col)
            tf2 = features2.add_feature(f2, color=col, border=col)
            self.diagram.cross_track_links.append(CrossLink(tf1, tf2, col, col))

    def _compose_crosslinks(self):
        if not self.diagram or not self.fsets: return
        print 'Finding crosslinks between genes in the cluster.'
        num_cids = len(self.clusters)
        for ci1, c1 in enumerate(self.clusters):
            if ci1 >= num_cids - 1: break
            features1 = self.fsets[ci1]
            ci2 = ci1 + 1
            features2 = self.fsets[ci2]
            c2 = self.clusters[ci2]
            print '%s vs %s' % (c1.sequence.description, c2.sequence.description)
            work = self.Work()
            work.start_work(self._blast_feature, c1.CDS, None, c1.sequence, c2.sequence)
            work.assemble(self._compose_crosslink, features1, features2)
            work.wait()

    def draw_clusters(self):
        print 'Creating cluster diagram.'
        # create diagram
        self.diagram = GenomeDiagram.Diagram(self.name)
        # add tracks
        max_len = max(len(c) for c in self.clusters)
        normal_color = colors.grey if self.no_color else colors.black
        mark_color = colors.black if self.no_color else colors.red
        for cluster in self.clusters:
            col = mark_color if cluster.marked else normal_color
            track = self.diagram.new_track(1,
                                           name=cluster.label,
                                           greytrack=1, height=0.4,
                                           greytrack_fontcolor=col,
                                           greytrack_labels=1,
                                           greytrack_fontsize=self.name_size,
                                           scale=False,
                                           start=0, end=max_len)
            self.fsets.append(track.new_set())
        # add crosslink features
        if self.add_crosslinks:
            self._compose_crosslinks()
        # add CDS-es
        for ci, cluster in enumerate(self.clusters):
            gene_id = 0
            for f in cluster.CDS:
                known = False
                fname = 'NONE'
                if gene_id < len(cluster.genes):
                    fname = cluster.genes[gene_id]
                    fcolor = self.colors.get(fname, self.default_color)
                if fname == 'NONE':
                    fname = self._feature_name(f, default='')
                    fcolor = self.default_color
                else:
                    #decapitalize gene names if they are marked as proteins
                    # fname = fname[:1].lower() + fname[1:] if fname else ''
                    known = True
                if self.no_color: fcolor = colors.color2bw(fcolor)
                self.fsets[ci].add_feature(f, sigil="BIGARROW",
                                           color=fcolor,
                                           border=fcolor if self.no_border else colors.black,
                                           name=fname, label=True,
                                           label_position="middle",
                                           label_size=self.gene_size if known else self.unknown_gene_size,
                                           label_color=colors.black if known else colors.grey,
                                           label_angle=self.gene_angle)
                gene_id += 1
        self.diagram.draw(format="linear", pagesize=self.pagesize, fragments=1,
                          start=0, end=max_len)
        for ptype in ('PDF', 'EPS', 'SVG'):
            dianame = '%s.%s' % (self.name, ptype.lower())
            print 'Saving: %s' % dianame
            self.diagram.write(dianame, ptype)
        print 'Done.'


