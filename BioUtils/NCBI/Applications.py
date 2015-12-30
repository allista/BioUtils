# coding=utf-8
# Copyright (C) 2015 Allis Tauri <allista@gmail.com>
# 
# BioUtils is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# indicator_gddccontrol is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
Created on Dec 29, 2015

@author: Allis Tauri <allista@gmail.com>
'''

from __future__ import print_function

from Bio.Application import AbstractCommandline, _Switch, _Option, _Argument

__docformat__ = "restructuredtext en"

class FormatDBCommandline(AbstractCommandline):
    """Commandline object for the NCBI Blast+ formatdb program.
    """
    def __init__(self, cmd="formatdb", **kwargs):
        isstr = lambda x: isinstance(x, str)
        self.parameters = [
            _Switch(["--help", "help"],
                    "show brief help on version and usage"),
            _Option(["-t", "title"],
                    "Title for database file [String]  Optional",
                    checker_function=isstr,
                    equate=False),
            _Option(["-i", "input"],
                    "Input file(s) for formatting [File In]  Optional",
                    checker_function=isstr,
                    filename=True,
                    equate=False),
            _Option(["-l", "logfile"],
                    "Logfile name: [File Out]  Optional. Default = formatdb.log",
                    checker_function=isstr,
                    filename=True,
                    equate=False),
            _Option(["-p", "protein"],
                    "Type of file: T - protein, F - nucleotide [T/F]  Optional. Default = T",
                    checker_function=isstr,
                    equate=False),
            _Option(["-o", "options"],
                    "Parse options: T - True: Parse SeqId and create indexes. "
                    "F - False: Do not parse SeqId. Do not create indexes. [T/F]  Optional. Default = F",
                    checker_function=isstr,
                    equate=False),
            _Option(["-a", "ans_input"],
                    "Input file is database in ASN.1 format (otherwise FASTA is "
                    "expected). T - True, F - False. [T/F]  Optional. Default = F",
                    checker_function=isstr,
                    equate=False),
            _Option(["-b", "bin_ans"],
                    "ASN.1 database in binary mode. T - binary, F - text mode. "
                    "[T/F]  Optional. Default = F",
                    checker_function=isstr,
                    equate=False),
            _Option(["-e", "seq_input"],
                    "Input is a Seq-entry [T/F]  Optional. Default = F",
                    checker_function=isstr,
                    equate=False),
            _Option(["-n", "name"],
                    "Base name for BLAST files [String]  Optional",
                    checker_function=isstr,
                    equate=False),
            _Option(["-v", "volume"],
                    "Database volume size in millions of letters [Integer]  Optional. Default = 4000",
                    checker_function=isstr,
                    equate=False),
            _Option(["-s", "sparse"],
                    "Create indexes limited only to accessions - sparse [T/F]  Optional. Default = F",
                    checker_function=isstr,
                    equate=False),
            _Option(["-V", "verbose"],
                    "Verbose: check for non-unique string ids in the database [T/F]  Optional. Default = F",
                    checker_function=isstr,
                    equate=False),
            _Option(["-L", "alias"],
                    "Create an alias file with this name. Use the gifile arg (below) if set to calculate db size; "
                    "use the BLAST db specified with -i (above) [File Out]  Optional",
                    checker_function=isstr,
                    equate=False),
            _Option(["-F", "gifile"],
                    "Gifile (file containing list of gi's) [File In]  Optional",
                    checker_function=isstr,
                    filename=True,
                    equate=False),
            _Option(["-B", "gi_out"],
                    "Binary Gifile produced from the Gifile specified above [File Out]  Optional",
                    checker_function=isstr,
                    equate=False),
            _Option(["-T", "taxid"],
                    "Taxid file to set the taxonomy ids in ASN.1 deflines [File In]  Optional",
                    checker_function=isstr,
                    equate=False),
        ]
        super(FormatDBCommandline, self).__init__(cmd, **kwargs)
