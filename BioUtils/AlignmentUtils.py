
'''
Created on Jul 20, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import re

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Applications import FastTreeCommandline

from BioUtils.SeqUtils import mktmp_fasta, SeqLoader, copy_attrs, num_fasta_records
from BioUtils.Tools.Text import FilenameParser, issingleletter
from BioUtils.Tools.Multiprocessing import cpu_count
from BioUtils.Tools.Misc import safe_unlink, run_cline, mktmp_name

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
    
    def index(self, base_sid):
        return self._db.get(base_sid, None)
    
    @property
    def rows(self): return xrange(len(self))
    
    @property
    def cols(self): return xrange(self.get_alignment_length())
    
    @property
    def rcols(self): return xrange(self.get_alignment_length()-1,-1,-1)
    
    def remove(self, *sids):
        rows = set(self.rows)
        for base_sid in sids:
            i = self.index(base_sid)
            if base_sid is None: continue
            rows.remove(i)
        rows = sorted(rows)
        return AlignmentExt([self[i] for i in rows], self._alphabet, self.annotations)
    
    def trim(self, gap = '-'):
        '''Trim alignment termina containing gaps in columns 
        and columns containing only gaps'''
        #trim margins
        lmargin = 0
        for i in self.cols:
            lmargin = i
            if gap not in self[:,i]: break
        if lmargin == self.get_alignment_length()-1:
            return AlignmentExt([])
        rmargin = -1
        for i in self.rcols:
            rmargin = i+1
            if gap not in self[:,i]: break
        ali = self[:,lmargin:rmargin]
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
        slices.append(slice(s, cols[-1]+1))
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
                 'phylip':  re.compile(r'.*\.phy$'),
    }
    
    @classmethod
    def load(cls, filename, schema=None):
        try: 
            return [AlignmentExt.from_msa(msa) for msa in
                    AlignIO.parse(filename, cls.schema(filename, schema))]
        except Exception, e:
            print 'Unable to load alignments from: %s\n%s' % (filename, str(e))
            return None
        
    @classmethod
    def load_first(cls, filename, schema=None):
        alis =  cls.load(filename)
        if not alis: return None
        return alis[0]
        
    @classmethod
    def save(cls, alignments, filename, schema=None):
        try: 
            AlignIO.write(alignments, filename, cls.schema(filename, schema))
            return True
        except Exception, e:
            print 'Unable to save alignments to: %s\n%s' % (filename, str(e))
            return False
        
    @classmethod
    def mktmp(cls, alignments, schema='fasta'):
        outfile = mktmp_name('.%s' % schema)
        cls.save(alignments, outfile, schema)
    
    @classmethod
    def align(cls, seq_records, outfile=None):
        '''Align given sequences
        @param seq_records: a list of SeqRecords objects
        @param outfile: a filename for the output alignment or None
        @return: if the outfile is none, return an AlignmentExt object;
        otherwise return True on success. In both cases return None on error.'''
        if not outfile: 
            outfile = mktmp_name('.aln.fasta')
            remove_out = True
        else: remove_out = False
        msafile = mktmp_fasta(seq_records)
        args = dict(thread=cpu_count, input=msafile)
        if len(seq_records) < 10000: 
            args['auto'] = True
        else: 
            args['parttree'] = True
            args['partsize'] = 1000
        ali = None
        if run_cline(MafftCommandline(**args), stdout=outfile):
            if remove_out:
                ali = AlignmentExt.from_msa(AlignIO.read(outfile, 'fasta'))
            else: ali = True
        if remove_out: safe_unlink(outfile)
        safe_unlink(msafile)
        return ali
    #end def
    
    @classmethod
    def align_translated(cls, seq_records, trans_table=1, force_trans_table=False):
        '''Translate given sequences, align them and untranslate.
        Return an alignment object'''
        trans_seqs = cls._translate(seq_records, trans_table, force_trans_table)
        alignment  = cls.align(trans_seqs)
        untrans_aligned_seqs = cls._untranslate(alignment, seq_records)
        return AlignmentExt(untrans_aligned_seqs)
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