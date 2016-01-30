'''
Created on Jan 30, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os
import time

from ..Taxonomy import Lineage

RETRIES = 3
PAUSE = 5

def KEGG_cmd(cmd, args, outfile, emsg, force=False):
    if not force and os.path.isfile(outfile): return True
    for _try in xrange(RETRIES):
        try: 
            request = cmd(args)
            with open(outfile, 'w') as out:
                out.write(request.read())
            return True
        except Exception, ex:
            print ex
            if _try < RETRIES-1: time.sleep(PAUSE)
            else: print emsg 
    return False
    

class KGene(object):
    '''Represents KEGG's gene database entry in the form organism:gen_id'''
    def __init__(self, org, gid, desc=''):
        self.org = org
        self.id  = gid
        self.description = desc
        self.request = "%s:%s" % (self.org, self.id) 
        self.filename = "%s-%s.kegg" % (self.org, self.id) 
    
    def __str__(self): return self.request
#end class

class Organism(object):
    '''Holds a KEGG's organism record as returned by /list/organism REST command'''
    def __init__(self, line=None):
        self.T = ""
        self.org = ""
        self.description = ""
        self.lineage = None
        if line is not None: self.parse(line)
        
    def __str__(self):
        return '%s: %s [%s]' % (self.org, self.description, self.lineage)
    
    def parse(self, line):
        fields = line.split('\t')
        if len(fields) < 4:
            raise ValueError("Organism line should have three fields divided by tabs")
        self.T = fields[0]
        self.org = fields[1]
        self.description = fields[2]
        self.lineage = Lineage(fields[3])
        
    def belongs2(self, lineage):
        return lineage.includes(self.lineage)
    
    def fix_lineage(self):
        genus = self.description.split()[0].lower()
        if genus in Lineage.known_taxons:
            if genus not in self.lineage: 
                self.lineage = Lineage.from_iter(list(self.lineage)+[genus])
            else:
                self.lineage = Lineage.from_iter(self.lineage[:self.lineage.index(genus)+1])
    
    @classmethod
    def from_record(cls, rec):
        org = cls()
        org.T = rec.id
        org.org = rec.id
        org.description = rec.annotations.get('organism', rec.name)
        org.lineage = Lineage.from_iter(rec.annotations.get('taxonomy'))
        return org
#end class

class Organisms(dict):
    '''Representation and parser for KEGG's organism list'''
    def __init__(self, filename=None):
        dict.__init__(self)
        self.parse(filename)
        
    def parse(self, filename):
        if not filename or not os.path.isfile(filename): return False
        print 'Parsing the list of organisms.'
        with open(filename) as inp:
            for line in inp:
                try: org = Organism(line)
                except ValueError, er:
                    print 'Unable to parse organism info: %s' % line
                    print er
                    continue
                self[org.org] = org
                if org.description: self[org.description] = org
        return True
    
    def fix_lineages(self):
        for org in self.values():
            org.fix_lineage()
#end class

class KEGGRecord(dict):
    '''Generic representation of the the KEGG file format'''
    def __init__(self):
        dict.__init__(self)
        self.entry = ""
        self.description = ""
        
    def __str__(self):
        return 'ENTRY: %s\t%s\n%s' % (self.entry, self.description, dict.__str__(self))
    
    @classmethod
    def parse(cls, filename):
        if not filename or not os.path.isfile(filename): return
        with open(filename) as inp:
            record = cls()
            keyword = ""
            for line in inp:
                if line[:3] == '///':
                    yield record
                    record = cls()
                    continue
                if line[:12] != "            ":
                    keyword = line[:12].strip()
                data = line[12:].strip()
                if keyword == "ENTRY":
                    words = data.split()
                    record.entry = words[1]
                    record.description = " ".join(words[1:])
                field = record.get(keyword, None)
                if field is None:
                    field = []
                    record[keyword] = field
                field.append(data)
#end class