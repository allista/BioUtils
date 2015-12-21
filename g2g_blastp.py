#!/usr/bin/python
# coding=utf-8
#
# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
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
Created on Dec 19, 2015

@author: Allis Tauri <allista@gmail.com>
'''

import os, re
from BioUtils.CommonTools import isatty
from BioUtils.Blast import LocalBlast

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC, NucleotideAlphabet, ProteinAlphabet

from DegenPrimer.MultiprocessingBase import MultiprocessingBase, cpu_count


if __name__ == '__main__':
    pass