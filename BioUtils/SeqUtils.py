# coding=utf-8

'''
Created on Dec 24, 2015

@author: Allis Tauri <allista@gmail.com>
'''

import os, sys, re
import tempfile, itertools
from BioUtils.Tools.Output import user_message
from copy import deepcopy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_alphabet, generic_dna, generic_protein, ProteinAlphabet, DNAAlphabet, NucleotideAlphabet
from Bio.Data.IUPACData import ambiguous_dna_letters, ambiguous_rna_letters, extended_protein_letters

from .Tools.Multiprocessing import MultiprocessingBase, ordered_shelved_results_assembler
from .Tools.tmpStorage import shelf_result, register_tmp_file, cleanup_file
from .Tools.Text import FilenameParser, random_text
from .Tools.Misc import mktmp_name, safe_unlink

re_type = type(re.compile(''))
isatty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()

class SeqLoader(MultiprocessingBase, FilenameParser):
    fasta_re = re.compile(r'.*\.f(asta|as|a|aa|fn|rn|na)$')
    
    schemas = {'fasta': fasta_re,
               'gb':    re.compile(r'.*\.gb(ff|k)$'),
               'embl':  re.compile(r'.*\.embl$'),
    }
    
    def __init__(self, abort_event):
        super(SeqLoader, self).__init__(abort_event)

    dna_letters = set(ambiguous_dna_letters.upper())
    rna_letters = set(ambiguous_rna_letters.upper())
    protein_letters = set(extended_protein_letters.upper())

    @classmethod
    def guess_alphabet(cls, seq):
        letters = set(seq[:10].upper()) # use just first 10 letters
        nletters = len(letters)
        if len(letters.intersection(cls.dna_letters)) == nletters:
            return generic_dna
        if len(letters.intersection(cls.protein_letters)) == nletters:
            return generic_protein
        return generic_alphabet

    @classmethod
    def correct_alphabet(cls, seq):
        cls._set_alphabet(seq, cls.guess_alphabet(seq))
        return seq

    @classmethod
    def _set_alphabet(cls, rec, abc):
        rec.seq.alphabet = abc
        return rec

    @classmethod
    def load_file(cls, filename, schema=None, guess_alphabet=False):
        if not os.path.isfile(filename):
            print 'No such file: %s' % filename 
            return None
        if not schema: schema = cls.guess_schema(filename) 
        try:
            if guess_alphabet:
                recs = list(cls.correct_alphabet(rec) for rec in SeqIO.parse(filename, schema))
            else:
                recs = list(SeqIO.parse(filename, schema))
            return recs
        except Exception, e:
            print 'Unable to parse %s as %s\n%s' % (filename, schema, str(e))
            return None
        
    def load_dir(self, dirname, schema=None, namefilter=None, flatten=True, guess_alphabet=False):
        if isinstance(namefilter, str):
            namefilter = re.compile(namefilter)
        if isinstance(namefilter, re_type):
            flt = namefilter.match
        elif hasattr(namefilter, '__call__'):
            flt = namefilter
        else: flt = lambda n: True
        files = [f for f in (os.path.join(dirname, fn) 
                             for fn in os.listdir(dirname) if flt(fn))
                 if os.path.isfile(f)]
        if not files:
            print 'No files found.'
            return None
        return self.load_files(files, schema, flatten, guess_alphabet)
    
    def load_files(self, files, schema=None, flatten=True, guess_alphabet=False):
        results = self.parallelize_work(1, self.load_file, files, schema, guess_alphabet)
        if not results: return None
        results = filter(lambda x: bool(x), results)
        if flatten: return list(itertools.chain.from_iterable(results))
        else: return results

class SeqView(object):
    '''Uses Bio.SeqIO.index_db with a temporary database to create a mixed
    dict/list picklable read-only SeqRecord container.
    '''
    
    def __init__(self, upper=False):
        self.upper  = upper
        self.dbname = None
        self.db     = None
        self._ids   = []
        self.tmp_db = False
        self.master = False
        
    def __str__(self):
        return '\n'.join(('SeqView:', 
                          'dbname:  %s' % self.dbname,
                          'tmp_db:  %s' % self.tmp_db,
                          'master:  %s' % self.master,
                          'upper:   %s' % self.upper,
                          'records: %d' % len(self._ids)))

    def load_dir(self, dirname, dbname=None, namefilter=None):
        if isinstance(namefilter, str):
            namefilter = re.compile(namefilter)
        if isinstance(namefilter, re_type):
            flt = namefilter.match
        elif hasattr(namefilter, '__call__'):
            flt = namefilter
        else: flt = lambda n: True
        files = [f for f in (os.path.join(dirname, fn)
                             for fn in os.listdir(dirname) if flt(fn))
                 if os.path.isfile(f)]
        if not files:
            print 'No files found in: %s' % dirname
            return False
        self.load(files, dbname)

    def load(self, files, dbname=None):
        if isinstance(files, basestring): files = [files]
        self.close()
        valid = []
        schemas = set()
        for filename in files:
            if isinstance(filename, bytes):
                filename = filename.decode('utf8')
            elif isinstance(filename, str):
                filename = unicode(filename)
            if not os.path.isfile(filename):
                print 'No such file: %s' % filename 
                continue
            schema = SeqLoader.guess_schema(filename)
            if not schema:
                print 'Unable to guess schema from filename: %s' % filename
                continue
            schemas.add(schema)
            valid.append(filename)
        if len(schemas) != 1:
            raise ValueError('All files should be of the same type, but %d types found: %s' % (len(schemas), ', '.join(schemas)))
        if not valid:
            print 'No valid files provided.'
            return False
        if not dbname:
            self.dbname = mktmp_name('_SeqView.db')
            safe_unlink(self.dbname)
            self.tmp_db = True
        else: self.dbname = dbname
        self.db = SeqIO.index_db(self.dbname, valid, schemas.pop())
        self._ids = tuple(sorted(self.db.keys()))
        self.master = True
        return bool(self)

    @classmethod
    def safe_load(cls, files):
        with user_message("Loading sequences..."):
            try:
                view = cls()
                view.load(files)
            except Exception as e:
                print str(e)
                return None
        if len(view) == 0:
            print 'No sequences were loaded from:\n%s' % files
            view.close()
            return None
        return view

    def reload(self, dbname):
        self.close()
        self.dbname = dbname
        self.db = SeqIO.index_db(self.dbname)
        self._ids = sorted(self.db.keys())
        self.tmp_db = False
        self.master = True
        
    def close(self):
        if self.master: 
            self.db.close()
            self.master = False
        self.db = None
        self._ids = []
        if self.dbname and self.tmp_db:
            cleanup_file(self.dbname)
            self.dbname = None
            self.tmp_db = False
            
    def __del__(self): self.close()
            
    def __nonzero__(self): return bool(self._ids)
        
    def __len__(self): return len(self._ids)
    
    def __iter__(self): return (self.db[key] for key in self._ids)
    
    def __getitem__(self, key):
        rec = None
        if isinstance(key, int):
            rec = self[self._ids[key]]
        elif isinstance(key, str):
            rec = self.db[key]
        elif isinstance(key, slice):
            rec = self.subview(self._ids[key])
        if rec is None: raise KeyError('Sequence %s not found' % key)
        return rec.upper() if self.upper else rec
        
    def get(self, key, default=None):
        try: return self[key]
        except: return default
    
    def iterkeys(self): return self._ids.__iter__()
    def keys(self): return self._ids
    def key(self, index): return self._ids[index]
    
    def itervalues(self, keys=None):
        if not keys: return self.__iter__()
        return (self[k] for k in keys)
        
    def values(self, keys=None):
        return list(self.itervalues(keys))
    
    def iteritems(self):
        for k in self.iterkeys(): yield (k, self[k])
        
    def items(self): return dict(self.iteritems())
        
    def subview(self, keys):
        v = SeqView()
        v.upper = self.upper
        v.dbname = self.dbname
        v.db = self.db
        v.tmp_db = False
        v._ids = keys
        return v
    
    def clone(self):
        return self.__deepcopy__(dict())
    
    def __deepcopy__(self, memo):
        return _unpickle_SeqView(self.dbname, deepcopy(self._ids, memo), self.upper)
    
    def __reduce__(self):
        return _unpickle_SeqView, (self.dbname, self._ids, self.upper)

def _unpickle_SeqView(dbname, ids, upper):
    v = SeqView(upper)
    v.dbname = dbname
    v.db = SeqIO.index_db(dbname)
    v.master = True
    v._ids = ids
    return v

    
class Translator(MultiprocessingBase):
    def __init__(self, abort_event):
        super(Translator, self).__init__(abort_event)
        
    @staticmethod
    @MultiprocessingBase.data_mapper
    @shelf_result
    def _translate_feature(fi, rec, table):
        f = rec.features[fi]
        srec = f.extract(rec)
        try: 
            tsec = srec.seq.translate(table)
            if tsec[-1] == '*': tsec = tsec[:-1]
            try: fid = f.qualifiers['locus_tag'][0]
            except KeyError: fid = '%s_f%d' % (rec.id, fi)
            trec = SeqRecord(tsec,
                             id=fid, name=rec.name,
                             description=rec.description,
                             annotations=rec.annotations)
            pf = SeqFeature(FeatureLocation(0, len(trec)), 
                            id=f.id, type='CDS',
                            qualifiers=f.qualifiers)
            trec.features.append(pf)
        except Exception, e:
            print e
            raise RuntimeError('Unable to translate: %s' % str(srec.seq))
        return trec

    @staticmethod
    def _translate_sixframe_seq(i, recs, table):
        rec = recs[i]
        try:
            tseq = rec.seq.translate(table)
            if tseq[-1] == '*': tseq = tseq[:-1]
            trec = SeqRecord(tseq,
                             id='FRAME%d'%i, name=rec.name,
                             description=rec.description,
                             annotations=rec.annotations)
        except Exception, e:
            print e
            raise RuntimeError('Unable to translate frame %d of: %s %s' % (i, rec.id, rec.description))
        return trec

    @staticmethod
    def translate(record, table):
        rem = len(record) % 3
        if rem > 0: record = record[:-rem]
        tseq = record.seq.translate(table)
        if tseq[-1] == '*': tseq = tseq[:-1]
        trec = SeqRecord(tseq,
                         id=record.id, name=record.name,
                         description=record.description,
                         annotations=record.annotations)
        return trec

    def translate_records(self, records, table):
        @MultiprocessingBase.data_mapper
        @shelf_result
        def worker(i, recs):
            return Translator.translate(recs[i], table)
        translation = [None]*len(records)
        work = self.Work()
        work.start_work(worker, range(len(records)), None, records)
        work.assemble(ordered_shelved_results_assembler, translation)
        if not work.wait(): return None
        return translation

    def translate_features(self, rec, features=None, table='Standard', join=False, gap='X' * 20):
        if features is None: features = rec.features
        translation = [None]*len(features)
        work = self.Work()
        work.start_work(self._translate_feature, features, None, rec, table)
        work.assemble(ordered_shelved_results_assembler, translation)
        if not work.wait(): return None
        if join: return cat_records(translation, gap=gap)
        return translation

    def _get_frame(self, rec, frame, revcomp=False):
        frec = rec[frame:]
        rem = len(frec) % 3
        if rem > 0: frec = frec[:-rem]
        frec.annotations = dict(rec.annotations)
        frec.dbxrefs = list(rec.dbxrefs)
        frec.annotations['start'] = frame
        frec.annotations['strand'] = -1 if revcomp else 1
        return frec

    def _get_six_frames(self, rec):
        if len(rec) < 6:
            raise ValueError('Sequence %s is too short for 6-frame translation. '
                             'Should be at least 6 letters long.' % rec.id)
        #prepare six frames
        revcomp = rec.reverse_complement()
        recs = [self._get_frame(rec, 0),
                self._get_frame(rec, 1),
                self._get_frame(rec, 2),
                self._get_frame(revcomp, 0, True),
                self._get_frame(revcomp, 1, True),
                self._get_frame(revcomp, 2, True)]
        return recs

    _frames = (0,1,2,3,4,5)

    def translate_six_frames(self, rec, table='Standard'):
        @MultiprocessingBase.data_mapper
        @shelf_result
        def worker(i, frames):
            return Translator._translate_sixframe_seq(i, frames, table)
        translation = [None] * 6
        work = self.Work()
        work.start_work(worker, self._frames, None, self._get_six_frames(rec))
        work.assemble(ordered_shelved_results_assembler, translation)
        if not work.wait(): return None
        return translation

    def translate_six_frames_single(self, rec, table='Standard'):
        '''Single-process variant of translate_six_frames'''
        recs = self._get_six_frames(rec)
        return [self._translate_sixframe_seq(i, recs, table) for i in self._frames]

    
def load_dir(abort_event, dirname, schema=None, namefilter=None, flatten=True, guess_alphabet=False):
    loader = SeqLoader(abort_event)
    return loader.load_dir(dirname, schema, namefilter, flatten, guess_alphabet)


def load_files(abort_event, filenames, schema=None, flatten=True, guess_alphabet=False):
    loader = SeqLoader(abort_event)
    return loader.load_files(filenames, schema, flatten, guess_alphabet)


def common_alphabet(seqs):
    alphabets = set()
    last = None
    for s in seqs:
        if isinstance(s.seq.alphabet, ProteinAlphabet):
            last = ProteinAlphabet
            alphabets.add(1)
        elif isinstance(s.seq.alphabet, NucleotideAlphabet):
            last = NucleotideAlphabet
            alphabets.add(2)
        else:
            last = type(s.seq.alphabet)
            alphabets.add(type)
        if len(alphabets) > 1:
            return False, None
    return True, last


def mktmp_fasta(rec, register=True):
    fd, fn = tempfile.mkstemp('.fas', 'fa')
    f = os.fdopen(fd, 'wb')
    SeqIO.write(rec, f, 'fasta')
    if register: register_tmp_file(fn)
    f.close()
    return fn
#end def

def num_fasta_records(filename):
    num_records = 0
    with open(filename) as inp:
        for l in inp:
            if l[0] == '>':
                num_records += 1
    return num_records

def cat_records(records, gap='X'*20):
    cat = None
    for t in records:
        if t: 
            if cat: cat += gap+t
            else: cat = t
    return cat
#end def


def unique_records(*records):
    ids = set()
    for recs in records:
        for r in recs:
            if r.id in ids: continue
            ids.add(r.id)
            yield  r


def safe_parse(filename, schema=None, alphabet=None):
    try: return SeqIO.parse(filename, SeqLoader.schema(filename, schema), alphabet)
    except Exception, e:
        print 'Unable to parse %s\n%s\n' % (filename, str(e))
        return []


def safe_write(records, filename, schema=None):
    try: SeqIO.write(records, filename, SeqLoader.schema(filename, schema))
    except Exception, e:
        print 'Unable to write records to %s\n%s\n' % (filename, str(e))
        return False
    return True

def backup_write(records, filename, schema=None):
    if os.path.isfile(filename):
        os.rename(filename, filename + '.back')
    return safe_write(records, filename)

def features_of_type(record, ftype):
    return [f for f in record.features if f.type == ftype]

def all_CDS(record):
    return [f for f in record.features if f.type == 'CDS']

def all_genes(record):
    return [f for f in record.features if f.type == 'gene']

def features_of_type_indexes(record, ftype):
    return [i for i, f in enumerate(record.features) if f.type == ftype]


def get_indexes_of_all_genes(rec):
    features = features_of_type_indexes(rec, 'CDS')
    if not features:
        features = features_of_type_indexes(rec, 'gene')
        if not features:
            print 'No gene/CDS features found in:\n%s %s' % (rec.id, rec.description)
            return None
    return features


def get_feature_indexes_by_qualifier(rec, qualifier, values):
    if not isinstance(values, (list, tuple)):
        values = [values]
    indexes = dict.fromkeys(values)
    for i, f in enumerate(rec.features):
        try: tag = f.qualifiers[qualifier][0]
        except: continue
        try: index = indexes[tag]
        except: continue
        if index is not None: continue
        indexes[tag] = i
    return [i for i in (indexes.get(tag) for tag in values) if i is not None]


def get_indexes_of_genes(rec, genes):
    return get_feature_indexes_by_qualifier(rec, 'locus_tag', genes)


def simple_rec(seq, sid, name='', description='', feature='', alphabet=generic_alphabet):
    rec = SeqRecord(Seq(seq, alphabet), id=sid, name=name, description=description)
    if feature: rec.features.append(simple_feature(0, len(rec), feature, feature))
    return rec


def simple_feature(start, end, fid='<unknown id>', ftype='misc_feature', qualifiers=None):
    '''Make a SeqFeature starting at start and ending at end.
    Indexes are zero-based, end is not included and if start > end
    the feature is constructed on the reverse-complement strand.'''
    if start > end: loc = FeatureLocation(end, start, -1)
    else: loc = FeatureLocation(start, end, 1)
    return SeqFeature(loc, type=ftype, id=fid, qualifiers=qualifiers)


_illegal = re.compile(r'[=/]')

def pretty_rec_name(rec):
#    print 'source: %s, id %s, desc %s' % (rec.annotations.get('source'), rec.id, rec.description)#test
    if rec.description is None: rec.description = ''
    name = (rec.annotations.get('source')
            or rec.description.replace(rec.id, '').strip() 
            or rec.id)
    name = _illegal.sub('-', name)
    return name

def rec_name_with_id(rec):
    return '%s (%s)' % (pretty_rec_name(rec), rec.id)

def filename_for_record(rec, ext='gb'):
    return '%s_%s.%s' % (pretty_rec_name(rec).replace(' ', '_').strip('.'), rec.id, ext)

def copy_attrs(frec, trec):
    trec.id = frec.id
    trec.name = frec.name
    trec.description = frec.description
    return trec


def random_DNA(length): return random_text(length, 'ATGC')
def random_DNA_rec(length, sid, name='', description='', feature=''):
    return simple_rec(random_DNA(length), sid, name, description, feature, alphabet=generic_dna)
    
