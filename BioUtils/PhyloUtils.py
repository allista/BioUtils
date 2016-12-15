"""
Created on Mar 19, 2016

@author: Allis Tauri <allista@gmail.com>
"""

import re
import os 
import dendropy as dp
from xml.dom import minidom

import math
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Applications import FastTreeCommandline

from .AlignmentUtils import AlignmentUtils
from .SeqUtils import num_fasta_records
from .Tools.Text import FilenameParser
from .Tools.Multiprocessing import MultiprocessingBase
from .Tools.Misc import run_cline, mktmp_name
from .Tools.Output import user_message
from .Tools.Debug import estr
from .Taxonomy import Lineage


class PhyloUtils(MultiprocessingBase, FilenameParser):
    schemas   = {'newick': re.compile(r'.*\.(tre|nwk|txt)$'),
    }

    reroot_midpoint = 'midpoint'
    reroot_unroot = 'unroot'
    reroot_specials = (reroot_midpoint, reroot_unroot)

    _organism_re = re.compile(r'\s\[.*?\]\.?$')
    
    @classmethod
    def build_fast_tree(cls, alignment, outfile, **kwargs):
        # type: (MultipleSeqAlignment, str) -> bool
        """Build an approximate-ML tree with fasttree.
        :param outfile: output filename
        :param alignment: alignment object or a filename of an alignment in fasta format"""
        if isinstance(alignment, MultipleSeqAlignment):
            alnfile = mktmp_name('.aln.fasta')
            AlignmentUtils.save(alignment, alnfile)
        elif isinstance(alignment, basestring):
            alnfile = alignment
        else: raise TypeError('Unsupported type of alignment argument: %s' % type(alignment))
        args = dict(input=alnfile, out=outfile, pseudo=1)
        if num_fasta_records(alnfile) >= 10000:
            args['fastest'] = True
            args['boot'] = 100
        args.update(kwargs)
        if not run_cline(FastTreeCommandline(**args), name=cls.strip_ext(outfile)):
            return False 
        print 'Done\n'
        return True
    
    @classmethod
    def _add_format(cls, annotatable, value, datatype='xsd:string', replace = False):
        if annotatable.annotations:
            for a in annotatable.annotations:
                if a.name == 'format':
                    if a.value and not replace:
                        a.value = ' %s %s;' % (a.value.strip(' ;'), value)
                    else:
                        a.value = ' %s;' % value
                    return
        annotatable.annotations.add_new('format', ' %s;' % value, datatype_hint=datatype)
        
    @classmethod
    def _collapse_node(cls, node, name, collapse_hard):
        leafs = node.leaf_nodes()
        if any(hasattr(l, 'dont_collapse') for l in leafs):
            return
        num_leafs = len(leafs)
        print 'Collapsing %d %s' % (num_leafs, str(name))
        if collapse_hard:
            node_size = int(math.ceil(2*math.sqrt(num_leafs*25/3)))
            node_hw = (node_size, node_size)
            for child in node.child_nodes():
                node.remove_child(child)
            node.label = '%s (%d collapsed)' % (name.capitalize(), num_leafs)
            cls._add_format(node, 'nh=%d nw=%d sh=2' % node_hw)
        else:
            for cn in node.preorder_internal_node_iter():
                cn.annotations.add_new('collapse', '', datatype_hint='xsd:string')

    @classmethod
    def _set_node_taxonomy(cls, node, top_lineage, parent_lineage,
                           collapse_taxa, collapse_at_level,
                           collapse_min_nodes, collapse_hard,
                           hide_support, hide_taxonomy,
                           lineage_colors):
        lineages = [l.edge.lineage for l in node.leaf_iter()
                    if hasattr(l.edge, 'lineage') and l.edge.lineage]
        lineage = Lineage.common(lineages)
        nlins = len(lineages)
        if lineage: lineage = lineage - top_lineage
        if lineage and nlins  > 1 and lineage.first in collapse_taxa:
            cls._collapse_node(node, lineage.first, collapse_hard)
        elif (collapse_at_level and
              nlins >= collapse_min_nodes and
              Lineage.samefirst(lineages, collapse_at_level)):
            cls._collapse_node(node, lineage[-1].capitalize(), collapse_hard)
        #annotate edges
        if lineage and lineage != parent_lineage:
            if node.parent_node:
                l = (lineage-parent_lineage)
                if l.last == Lineage.unknown: l = l[:-1]
                if l:
                    if node.is_internal():
                        node.edge.label = l.capitalize()
                        if hide_support:
                            cls._add_format(node, 'lv=0') #FIXME: does not work
                        if hide_taxonomy:
                            cls._add_format(node.edge, 'lv=0')
                    else:
                        if l.first and node.label.startswith(l.first.capitalize()): l = l[1:]
                        if l.last and node.label.startswith(l.last.capitalize()): l = l[:-1]
                        if l and not hide_taxonomy: node.label += ' [%s]' % l.capitalize()
            #colorize edges
            if lineage_colors:
                col = None
                for l in lineage[::-1]:
                    col = lineage_colors.get(l)
                    if col: break
                #construct color meta for Dendroscope
                if col:
                    value = 'fg=%d %d %d' % col
                    if hide_taxonomy: value += ' lv=0'
                    for cn in node.preorder_iter():
                        cls._add_format(cn.edge, value, replace=True)
        for child in node.child_node_iter():
            cls._set_node_taxonomy(child, top_lineage, lineage,
                                   collapse_taxa, collapse_at_level,
                                   collapse_min_nodes, collapse_hard,
                                   hide_support, hide_taxonomy,
                                   lineage_colors)

    @classmethod
    def parse_colors(cls, colors_string):
        # parse colors
        try: from reportlab.lib import colors
        except ImportError:
            raise NotImplementedError('ReportLab is not installed. Unable to parse html colors.')
        parsed = dict()
        for col in colors_string.split(';'):
            lc = col.strip().split(':')
            if len(lc) == 2:
                col = lc[1]
                if col.startswith('#'):
                    col = colors.HexColor(col)
                else: col = getattr(colors, col, None)
                if col: parsed[lc[0]] = col.bitmap_rgb()
        return parsed

    @classmethod
    def load(cls, treefile, schema=None):
        if not os.path.isfile(treefile):
            print 'No tree file found.'
            return None
        try: return dp.Tree.get(path=treefile, schema=cls.schema(treefile, schema))
        except Exception, e:
            print 'Error while loading %s:\n%s' % (treefile, estr(e))
            return None
    
    @classmethod
    def annotate_tree(cls, treefile, organisms, outfile=None, schema=None, **kwargs):
        # type: (str, Taxonomy.Organisms, str, str, dict) -> bool
        '''
        Annotate input tree with taxonomy information using edge labels and colors.
        :param treefile : a file containing the tree to be annotated
        :param organisms: organisms database
        :param outfile: optional : basename for output file; Note: the last extension will be stripped
        :param schema: data format of the treefile

        Accepted kwargs:
        :param beautify leafs: bool (True) : replaces IDs in leafs' labels with organism names
        :param mark_leafs: list(str) : mark nodes with the specified labels in bold
        :param collapse_taxa : list(str) : collapses subtrees belonging to given taxa
        :param collapse_at_level : int (0) : collapses subtrees of the taxonomy of this level (1 - prokaryotes, 2 - bacteria...)
        :param collapse_hard: bool (False) : removes collapsed subtrees, leaving a single node; otherwise collapsed
        subtrees are displayed in Dendroscope as trapezium nodes
        :param collapse_min_nodes : int (3) : only collapse subtrees with number of leafs greater or equal than this
        :param min_support : float (0: disabled) : nodes with support less that this will be removed from the tree, children being relinked to parents
        :param reroot_at : string ('') : reroot the tree at specified leaf; special value 'midpoint' reroots the at midpoint; special value 'unroot' unroots the tree
        :param lineage_colors : dict : a dictionary of colors as (r, g, b) tuples with lowercase taxons as keys; special value 'auto' causes to automatically assign colors
        :param top_lineage: a Lineage object to be subtracted from lineages of organisms on the tree; if not provided, it is computed automatically
        :param hide_support : bool (False) : if True, internal node labels with support values are hidden in Dendroscope
        :param hide_taxonomy : bool (False) : if True, edge labels with taxonomy are hidden in Dendroscope
        '''
        with user_message('Processing tree file...', '\n'):
            tree = cls.load(treefile, schema)
            if not tree:
                print 'No tree loaded.'
                return False
            #need to get the root before beautifying
            new_root = None
            root_name = kwargs.pop('reroot_at', '')
            min_support = kwargs.pop('min_support', False)
            beautify_leafs = kwargs.pop('beautify_leafs', False)
            mark_leafs = kwargs.pop('mark_leafs', set())
            for leaf in tree.leaf_node_iter():
                label = leaf.taxon.label.replace(' ', '_')
                if label in mark_leafs:
                    cls._add_format(leaf, " x=0.0 y=0.0  ft='Ubuntu-BOLDITALIC-14' ll=7;")
                    setattr(leaf, "dont_collapse", True)
                if not new_root and label == root_name: new_root = leaf
                org = organisms.get(label)
                if not org: continue
                leaf.edge.lineage = org.lineage
                if beautify_leafs:
                    leaf.taxon.label = '%s' % (org.description or org.id)
                    if org.id != leaf.taxon.label:
                        leaf.taxon.label += ' (%s)' % org.id
                    leaf.taxon.label = leaf.taxon.label.replace('_', ' ')
                leaf.label = leaf.taxon.label
        #reroot the tree before traversing
        if root_name == cls.reroot_unroot:
            with user_message('Unrooting tree...'):
                tree.deroot()
        if root_name == cls.reroot_midpoint:
            with user_message('Rerooting tree at midpoint...'):
                tree.reroot_at_midpoint(update_bipartitions=True)
                tree.seed_node.label = 'midpoint root'
                cls._add_format(tree.seed_node, 'nh=10 nw=10 sh=4')
        elif new_root:
            with user_message('Rerooting tree at %s...' % root_name):
                tree.to_outgroup_position(new_root, update_bipartitions=True)
        else: print 'Node for rerooting not found: %s' % root_name
        #annotate the tree
        with user_message('Adding taxonomy information to the tree...', '\n'):
            top_lineage = kwargs.pop('top_lineage', None) 
            if not isinstance(top_lineage, Lineage): top_lineage = organisms.common_lineage
            colors = kwargs.pop('lineage_colors', None)
            if colors == 'auto': pass#TODO
            cls._set_node_taxonomy(tree.seed_node, top_lineage, None,
                                   kwargs.pop('collapse_taxa', []), 
                                   kwargs.pop('collapse_at_level', 0),
                                   kwargs.pop('collapse_min_nodes', 3), 
                                   kwargs.pop('collapse_hard', False),
                                   kwargs.pop('hide_support', False),
                                   kwargs.pop('hide_taxonomy', False),
                                   colors)
        if min_support:
            with user_message('Collapsing nodes with low support...'):
                for node in tree.postorder_internal_node_iter(exclude_seed_node=True):
                    try: support = float(node.label)
                    except(ValueError, TypeError): continue
                    if support < min_support and node.edge:
                        node.edge.collapse(adjust_collapsed_head_children_edge_lengths=True)
        with user_message('Saving resulting tree...'):
            if not outfile: outfile = cls.strip_ext(treefile)+'.out'
            xtreefile = outfile+'.nexml'
            tree.write(path=outfile+'.tre', schema='newick')
            tree.write(path=xtreefile, schema='nexml')
            with open(outfile+'.dot', 'w') as out:
                tree.write_as_dot(out, edge_formatter=lambda e: e.label or '')
        with user_message('Tuning nexml file for Dendroscope...'):
            cls._postprocess_nexml(xtreefile)
        return True
    
    @classmethod
    def replace_node_labels(cls, treefile, labels, schema=None, outfile=None):
        '''Reads a tree from file and replaces node labels
        according to provided mapping. The modified tree is
        returned as DendroPy.Tree object or is written to the 
        provided output file.
        :param labels: dict, replacement table
        :param outfile: the name of the file to write the modified tree
        '''
        with user_message('Loading tree file...', '\n'):
            tree = cls.load(treefile, schema)
            if not tree:
                print 'No tree loaded.'
                return None
        with user_message('Processing tree...', '\n'):
            for leaf in tree.leaf_node_iter():
                label = leaf.taxon.label.replace(' ', '_')
                if label in labels:
                    leaf.taxon.label = labels[label]
        if outfile: tree.write(path=outfile, schema=cls.schema(treefile, schema))
        return tree
    
    @staticmethod
    def _add_meta(e, content, datatype, prop, i):
        meta = minidom.Element('meta')
        attrs = {'content': content, 'datatype': 'xsd:%s'%datatype, 'id':'imeta%d'%i, 'property': prop, 'xsi:type':'nex:LiteralMeta'}
        for attr in attrs:  meta.setAttribute(attr, attrs[attr])
        if e.hasChildNodes(): e.insertBefore(meta, e.childNodes[0])
        else: e.appendChild(meta)
    
    @classmethod
    def _postprocess_nexml(cls, nexml_file):
        d = minidom.parse(nexml_file)
        #get collapsed nodes' ids
        collapsed = []
        for n in d.getElementsByTagName('meta'):
            lab = n.getAttribute('property')
            if not lab: continue
            if lab == 'dendropy:collapse' and n.parentNode:
                collapsed.append(n.parentNode.getAttribute('id'))
                n.parentNode.childNodes.remove(n)
            elif lab == 'dendropy:format':
                n.setAttribute('property', 'format') 
        #add the header
        for i, t in enumerate(d.getElementsByTagName('tree')):
            t.setAttribute('about', '#'+t.getAttribute('id'))
            t.setAttribute('label', '[%d]'%(i+1))
            t.setAttribute('xmlns:embellished', '')
            t.setAttribute('xmlns:collapsed', '')
            t.setAttribute('xmlns:drawer', '')
            t.setAttribute('xmlns:toscale', '')
            t.setAttribute('xmlns:radiallabels', '')
            t.setAttribute('xmlns:sparselabels', '')
            t.setAttribute('xmlns:root', '')
            cls._add_meta(t, "true", 'boolean', 'embellished', 2)
            cls._add_meta(t, "RectangularPhylogram", 'string', 'drawer', 1)
            cls._add_meta(t, "true", 'boolean', 'toscale', 2)
            cls._add_meta(t, "true", 'boolean', 'sparselabels', 3)
            cls._add_meta(t, "false", 'boolean', 'radiallabels', 4)
            cls._add_meta(t, " ".join(collapsed), 'string', 'collapsed', 5)
            cls._add_meta(t, " nh=2 nw=2 fg=0 0 0 bg=255 255 255 w=1 sh=0 fx=1 lc=0 0 0 lk=null ft='Ubuntu-ITALIC-13' lx=0 ly=0 ll=3 lv=1;", 
                           'string', 'default_node_format', 6)
            cls._add_meta(t, " fg=0 0 0 w=1 sh=1 dr=1 lc=0 0 0 lk=null ft='Ubuntu-ITALIC-13' lx=0 ly=0 ll=11 lv=1;", 
                           'string', 'default_edge_format', 7)
            cls._add_meta(t, t.getElementsByTagName('node')[0].getAttribute('id'), 'string', 'root', 8)
        tags = ('node', 'edge')
        for i, n in enumerate(d.getElementsByTagName('*')):
            if not n.tagName in tags: continue
            n.setAttribute('about', '#'+n.getAttribute('id'))
            no_format = True
            for sn in n.getElementsByTagName('meta'):
                if sn.getAttribute('property') == 'format':
                    no_format = False
                    break
            if no_format: cls._add_meta(n, ' ;', 'string', 'format', 100+i)
        with open(nexml_file, 'w') as out:
            out.write(d.toxml(encoding='utf-8'))
