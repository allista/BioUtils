
'''
Created on Jul 20, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import re

from StringIO import StringIO

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment

from BioUtils.SeqUtils import mktmp_fasta, SeqLoader, copy_attrs
from BioUtils.Tools.Text import FilenameParser, issingleletter
from BioUtils.Tools.Multiprocessing import cpu_count
from BioUtils.Tools.Misc import safe_unlink

class AlignmentExt(MultipleSeqAlignment):
    def __init__(self, records, alphabet=None, annotations=None):
        MultipleSeqAlignment.__init__(self, records, alphabet, annotations)
        self._db = dict((self[i].id, i) for i in xrange(len(self)))
        
    def __getitem__(self, index):
        if isinstance(index, str):
            i = self._db.get(index, None)
            if i is None: raise KeyError(index)
            return self[i]
        ret = MultipleSeqAlignment.__getitem__(self, index)
        if isinstance(ret, MultipleSeqAlignment):
            return self.from_msa(ret)
        return ret
    
    def __add__(self, other):
        return self.from_msa(MultipleSeqAlignment.__add__(self, other))
    
    def index(self, sid):
        return self._db.get(sid, None)
    
    @property
    def rows(self): return xrange(len(self))
    
    @property
    def cols(self): return xrange(self.get_alignment_length())
    
    @property
    def rcols(self): return xrange(self.get_alignment_length()-1,-1,-1)
    
    def remove(self, *sids):
        rows = set(self.rows)
        for sid in sids:
            i = self.index(sid)
            if sid is None: continue
            rows.remove(i)
        rows = sorted(rows)
        return AlignmentExt([self[i] for i in rows], self._alphabet, self.annotations)
    
    def trim(self, gap = '-'):
        '''Trim alignment termina containing gaps in columns 
        and columns containing only gaps'''
        #trim margins
        rmargin = 0
        for i in self.cols:
            rmargin = i
            if gap not in self[:,i]: break
        if rmargin == self.get_alignment_length()-1:
            return AlignmentExt([])
        lmargin = -1
        for i in self.rcols:
            lmargin = i
            if gap not in self[:,i]: break
        ali = self[:,rmargin:lmargin]
        #remove gaps-only
        cols = set(ali.cols)
        for i in ali.cols:
            if issingleletter(ali[:,i], gap):
                cols.remove(i)
        cols = sorted(cols)
        slices = []
        s = 0
        for i in xrange(len(cols)-1):
            if cols[i+1]-cols[i] > 1:
                slices.append(slice(s, cols[i]+1))
                s = cols[i+1]
        if not slices: 
            return AlignmentExt([])
        final = ali[:,slices[0]]
        for s in slices[1:]: final += ali[:,s]
        return final
    
    def reverse_complement(self):
        return AlignmentExt((copy_attrs(s, s.reverse_complement()) 
                             for s in self), self._alphabet)
    
    @classmethod
    def from_msa(cls, msa):
        return cls(msa, msa._alphabet, msa.annotations)

class AlignmentUtils(FilenameParser):
    schemas   = {'fasta':   SeqLoader.fasta_re,
                 'clustal': re.compile(r'.*\.aln$'),
                 'nexus':   re.compile(r'.*\.nex$'),
    }
    
    @classmethod
    def save(cls, alignments, filename, schema='clustal'):
        try: AlignIO.write(alignments, filename, cls.schema(filename, schema))
        except Exception, e:
            print 'Unable to save alignments to: %s' % filename
            print e
            
    @classmethod
    def load(cls, filename, schema=None):
        try: return list(AlignIO.parse(filename, cls.schema(filename, schema)))
        except Exception, e:
            print 'Unable to load alignments from: %s' % filename
            print e
    
    @classmethod
    def align(cls, seq_records, **kwargs):
        '''Align given sequences, return an alignment object'''
        try:
            msafile = mktmp_fasta(seq_records)
            args = dict(thread=cpu_count, input=msafile)
            if len(seq_records) < 10000: args['auto'] = True
            else: 
                args['parttree'] = True
                args['partsize'] = 1000
            cline = MafftCommandline(**args)
            print 'Running: %s' % cline
            out, err = cline()
            print err
            return AlignIO.read(StringIO(out), "fasta")
        except Exception, e:
            print 'AlignmentUtils.align: unable to align sequences:\n   %s' % str(e)
            return None
        finally: safe_unlink(msafile)
    #end def
    
    @classmethod
    def align_translated(cls, seq_records, trans_table=1, force_trans_table=False):
        '''Translate given sequences, align them and untranslate.
        Return an alignment object'''
        trans_seqs = cls._translate(seq_records, trans_table, force_trans_table)
        alignment  = cls.align(trans_seqs)
        untrans_aligned_seqs = cls._untranslate(alignment, seq_records)
        return MultipleSeqAlignment(untrans_aligned_seqs)
    #end def
    
    @classmethod
    def _translate(self, seq_records, trans_table=1, force_trans_table=False):
        '''Translate given sequences using trans_table. If some sequence feature 
        contains a translation table number, use it instead. If force_trans_table is True, 
        always use trans_table.'''
        trans_seqs =  list()
        for record in seq_records:
            #set the alphabet to DNA
            try:
                record.seq.alphabet = IUPAC.ambiguous_dna
            except Exception,e:
                #provided sequences SHOULD be DNA ones
                raise ValueError('AlignmentUtils.translate: unable to set alphabet of the sequence:\n%s\n%s' \
                                 % (record.format('fasta'), e.message))        
            #determine translation table
            translation_table = -1
            if force_trans_table:
                #force a translation table
                translation_table = trans_table
            else:
                #see if a translation table is defined in qualifiers
                for feature in record.features:
                    try:
                        translation_table = int(feature.qualifiers['transl_table'][0])
                        break
                    except: pass
                if translation_table < 0: translation_table = trans_table
            #do a translation
            trans_seq = record.seq.translate(table=translation_table, stop_symbol="X")
            trans_seq_rec = SeqRecord(trans_seq, id=record.id)
            trans_seq_rec.name = record.name
            trans_seq_rec.description = record.description
            trans_seqs.append(trans_seq_rec)
        return trans_seqs
    #end def    
    
    @classmethod
    def _untranslate(self, trans_seq_records, orig_seq_records):
        '''Untranslate given protein sequences using original nucleotide sequences, 
        thus conserving original triplets.'''
        #make a dictionary of sequences
        seq_dict = dict()
        for seq in orig_seq_records:
            seq_dict[seq.id] = seq
        #untranslate sequences
        untrans_seqs = list()
        for record in trans_seq_records:
            orig_seq = str(seq_dict[record.id].seq)
            rev_seq = ""
            l = 0
            for letter in record.seq:
                if letter == "-": rev_seq += "-"*3
                else:
                    if l+3 < len(orig_seq):
                        rev_seq += orig_seq[l:l+3]
                    else: rev_seq += orig_seq[l:]
                    l += 3
            if l < len(orig_seq): rev_seq += orig_seq[l:]
            untrans_seqs.append(SeqRecord(Seq(rev_seq, IUPAC.ambiguous_dna), id=record.id))
            #check results
            if rev_seq.replace("-", "") != orig_seq:
                raise Exception('AlignmentUtils.untranslate: resulting sequence differs from the original.')
        return untrans_seqs
    #end def
#end class