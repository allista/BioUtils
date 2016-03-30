#!/usr/bin/python
# coding=utf-8

'''
Created on Jul 9, 2012

@author: Allis Tauri <allista@gmail.com>
'''

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
import os
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

from distutils.core import setup
setup(name='BioUtils',
      version='1.4',
      description='A utility library and a set of scripts for everyday bioinformatic routine.',
      long_description=read('README'),
      license='GPL-3',
      author='Allis Tauri',
      author_email='allista@gmail.com',
      #url='https://launchpad.net/degenprimer',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Programming Language :: Python'],
#      install_require=['biopython>=1.6', 'dendropy>=4.0', 'numpy', 'scipy'],
      packages=['BioUtils', 'BioUtils.HMMER3', 'BioUtils.KEGG', 'BioUtils.NCBI', 'BioUtils.Tools'],
      scripts=['seqconv',
               'align_as_amino',
               'gaps_encoder',
               'parse_chrom',
               'split_fasta',
               'alien2gb',
               'gbfix',
               'rename_genomes',
               'analyze_ec',
               'split_consensus',
               'nameconv'
               ],
      #data_files=[('share/applications', ['DegenPrimer.desktop'])]
      )