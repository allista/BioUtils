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
Created on Aug 1, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from __future__ import print_function

from Bio.Application import AbstractCommandline, _Switch, _Option, _Argument

__docformat__ = "restructuredtext en"

class _HMMCommandlineBase(AbstractCommandline):
    """Base Commandline object for hmm commands from HMMER3.

    This is provided for subclassing, it deals with shared options
    common to all the hmm tools (hmmsearch, hmmscan, etc).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Core:
            _Switch(["-h", "h"],
                    "show brief help on version and usage"),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        super(_HMMCommandlineBase, self).__init__(cmd, **kwargs)


class _HMMCpuCommandlineBase(_HMMCommandlineBase):
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
        #Options controlling output:
        _Option(["--cpu", "cpu"],
                    "number of parallel CPU workers to use for multithreads",
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        super(_HMMCpuCommandlineBase, self).__init__(cmd, **kwargs)


class _HMMOutCommandlineBase(_HMMCommandlineBase):
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
        #Options controlling output:
        _Option(["-o", "o"],
                    "direct output to file <f>, not stdout",
                    filename=True,
                    equate=False),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        super(_HMMOutCommandlineBase, self).__init__(cmd, **kwargs)


class _HMMSeedCommandlineBase(_HMMOutCommandlineBase):
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
        #Other expert options:
        _Option(["--seed", "seed"],
                    "set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]",
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        super(_HMMSeedCommandlineBase, self).__init__(cmd, **kwargs)
        

class _HMMSearchCommandlineBase(_HMMSeedCommandlineBase, _HMMCpuCommandlineBase):
    """Base Commandline object for hmm search commands from HMMER3:
    hmmsearch, hmmscan
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
        #Options controlling output:
        _Option(["-tblout", "tblout"],
                    "save parseable table of per-sequence hits to file",
                    filename=True,
                    equate=False),
        _Option(["--domtblout", "domtblout"],
                    "save parseable table of per-domain hits to file",
                    filename=True,
                    equate=False),
        _Option(["--pfamtblout", "pfamtblout"],
                    "save table of hits and domains to file, in Pfam format",
                    filename=True,
                    equate=False),
        _Switch(["--acc", "acc"],
                "prefer accessions over names in output"),
        _Switch(["--noali", "noali"],
                "don't output alignments, so output is smaller"),
        _Switch(["--notextw", "notextw"],
                "unlimit ASCII text output line width"),
        _Option(["--textw", "textw"],
                    "set max width of ASCII text output lines  [120]  (n>=120)",
                    checker_function=lambda x: isinstance(x, int) and x >= 120,
                    equate=False),
        #Options controlling reporting thresholds:
        _Option(["-E", "E"],
                    "report sequences/profiles <= this E-value threshold in output  [10.0]  (x>0)",
                    checker_function=lambda x: isinstance(x, float) and x > 0,
                    equate=False),
        _Option(["-T", "T"],
                    "report sequences/profiles >= this score threshold in output",
                    equate=False),
        _Option(["--domE", "domE"],
                    "report domains <= this E-value threshold in output  [10.0]  (x>0)",
                    checker_function=lambda x: isinstance(x, float) and x > 0,
                    equate=False),
        _Option(["--domT", "domT"],
                    "report domains >= this score cutoff in output",
                    equate=False),
        #Options controlling model-specific thresholding:
        _Switch(["--cut_ga", "cut_ga"],
                "use profile's/sequence's GA gathering cutoffs to set all thresholding"),
        _Switch(["--cut_nc", "cut_tc"],
                "use profile's/sequence's NC gathering cutoffs to set all thresholding"),
        _Switch(["--cut_tc", "cut_nc"],
                "use profile's/sequence's TC gathering cutoffs to set all thresholding"),
        #Options controlling inclusion (significance) thresholds:
        _Option(["-incE", "incE"],
                    "consider sequences/profiles <= this E-value threshold as significant",
                    equate=False),
        _Option(["-incT", "incT"],
                    "consider sequences/profiles >= this score threshold as significant",
                    equate=False),
        _Option(["--incdomE", "incdomE"],
                    "consider domains <= this E-value threshold as significant",
                    equate=False),
        _Option(["--incdomT", "incdomT"],
                    "consider domains >= this score threshold as significant",
                    equate=False),
        #Options controlling acceleration heuristics:
        _Switch(["--max", "max"],
                "Turn all heuristic filters off (less speed, more power)"),
        _Option(["-F1", "F1"],
                    "Stage 1 (MSV) threshold: promote hits w/ P <= F1  [0.02]",
                    equate=False),
        _Option(["-F2", "F2"],
                    "Stage 2 (Vit) threshold: promote hits w/ P <= F2  [1e-3]",
                    equate=False),
        _Option(["--F3", "F3"],
                    "Stage 3 (Fwd) threshold: promote hits w/ P <= F3  [1e-5]",
                    equate=False),
        _Switch(["--nobias", "nobias"],
                "turn off composition bias filter"),
        #Other expert options:
        _Switch(["--nonull2", "nonull2"],
                "turn off biased composition score corrections"),
        _Option(["-Z", "Z"],
                    "set # of comparisons done, for E-value calculation",
                    equate=False),
        _Option(["-domZ", "domZ"],
                    "set # of significant seqs, for domain E-value calculation",
                    equate=False),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        super(_HMMSearchCommandlineBase, self).__init__(cmd, **kwargs)


class HMMSearchCommandline(_HMMSearchCommandlineBase):
    def __init__(self, cmd="hmmsearch", **kwargs):
        self.parameters = [
            #Options controlling output:
            _Option(["-A", "A"],
                    "save multiple alignment of all hits to file",
                    filename=True,
                    equate=False),
            #Other expert options:
            _Option(["-tformat", "tformat"],
                    "assert target <seqfile> is in format <s>: no autodetection",
                    checker_function=lambda s: isinstance(s, str),
                    equate=False),
            #Arguments
            _Argument(['hmmfile'],
                      "A path to a HMM profile",
                      filename=True, is_required=True),
            _Argument(['seqdb'],
                      "A path to a sequence database file to search in",
                      filename=True, is_required=True)
        ]
        super(HMMSearchCommandline, self).__init__(cmd, **kwargs)


class HMMScanCommandline(_HMMSearchCommandlineBase):
    def __init__(self, cmd="hmmscan", **kwargs):
        self.parameters = [
            #Other expert options:
            _Switch(["--daemon", "daemon"],
                    "run program as a daemon"),
            #Arguments
            _Argument(['hmmdb'],
                      "A path to a HMM profile database",
                      filename=True, is_required=True),
            _Argument(['seqfile'],
                      "A path to a sequence to search in",
                      filename=True, is_required=True)
        ]
        super(HMMScanCommandline, self).__init__(cmd, **kwargs)
        

class HMMBuildCommandline(_HMMSeedCommandlineBase, _HMMCpuCommandlineBase):
    def __init__(self, cmd="hmmbuild", **kwargs):
        self.parameters = [
        #Basic options:
        _Option(["-n", "n"],
                    "name the HMM",
                    equate=False),
        _Option(["--O", "O"],
                    "resave annotated, possibly modified MSA to file <f>",
                    filename=True,
                    equate=False),
        #Options for selecting alphabet rather than guessing it:         
        _Switch(["--amino", "amino"],
                "input alignment is protein sequence data"),
        _Switch(["--dna", "dna"],
                "input alignment is DNA sequence data"),
        _Switch(["--rna", "rna"],
                "input alignment is RNA sequence data"),
        #Alternative model construction strategies:
        _Switch(["--fast", "fast"],
                "assign cols w/ >= symfrac residues as consensus  [default]"),
        _Switch(["--hand", "hand"],
                "manual construction (requires reference annotation)"),
        _Option(["--symfrac", "symfrac"],
                    "sets sym fraction controlling --fast construction  [0.5]",
                    equate=False),
        _Option(["--fragthresh", "fragthresh"],
                    "if L <= x*alen, tag sequence as a fragment  [0.5]",
                    equate=False),
        #Alternative relative sequence weighting strategies:
        _Switch(["--wpb", "wpb"],
                "Henikoff position-based weights  [default]"),
        _Switch(["--wgsc", "wgsc"],
                "Gerstein/Sonnhammer/Chothia tree weights"),
        _Switch(["--wblosum", "wblosum"],
                "Henikoff simple filter weights"),
        _Switch(["--wnone", "wnone"],
                "don't do any relative weighting; set all to 1"),
        _Switch(["--wgiven", "wgiven"],
                "use weights as given in MSA file"),
        _Option(["--wid", "wid"],
                    "for --wblosum: set identity cutoff  [0.62]  (0<=x<=1)",
                    checker_function=lambda x: isinstance(x, float) and 0<=x<=1,
                    equate=False),
        #Alternative effective sequence weighting strategies:
        _Switch(["--eent", "eent"],
                "adjust eff seq # to achieve relative entropy target  [default]"),
        _Switch(["--eclust", "eclust"],
                "eff seq # is # of single linkage clusters"),
        _Switch(["--enone", "enone"],
                "no effective seq # weighting: just use nseq"),
        _Option(["-eset", "eset"],
                    "set eff seq # for all models to <x>",
                    equate=False),
        _Option(["-ere", "ere"],
                    "for --eent: set minimum rel entropy/position to <x>",
                    equate=False),
        _Option(["--esigma", "esigma"],
                    "for --eent: set sigma param to <x>  [45.0]",
                    equate=False),
        _Option(["--eid", "eid"],
                    "for --eclust: set fractional identity cutoff to <x>  [0.62]",
                    equate=False),
        #Alternative prior strategies:
        _Switch(["--pnone", "pnone"],
                "don't use any prior; parameters are frequencies"),
        _Switch(["--plaplace", "plaplace"],
                "use a Laplace +1 prior"),
        #Handling single sequence inputs:
        _Switch(["--singlemx", "singlemx"],
                "use substitution score matrix for single-sequence inputs"),
        _Option(["-popen", "popen"],
                    "gap open probability (with --singlemx)",
                    equate=False),
        _Option(["-pextend", "pextend"],
                    "gap extend probability (with --singlemx)",
                    equate=False),
        _Option(["-mx", "mx"],
                    "substitution score matrix (built-in matrices, with --singlemx)",
                    equate=False),
        _Option(["-mxfile", "mxfile"],
                    "read substitution score matrix from file <f> (with --singlemx)",
                    equate=False),
        #Control of E-value calibration:
        _Option(["-EmL", "EmL"],
                    "length of sequences for MSV Gumbel mu fit  [200]  (n>0)",
                    checker_function=lambda x: isinstance(x, int) and x>0,
                    equate=False),
        _Option(["-EmN", "EmN"],
                    "number of sequences for MSV Gumbel mu fit  [200]  (n>0)",
                    checker_function=lambda x: isinstance(x, int) and x>0,
                    equate=False),
        _Option(["-EvL", "EvL"],
                    "length of sequences for Viterbi Gumbel mu fit  [200]  (n>0)",
                    checker_function=lambda x: isinstance(x, int) and x>0,
                    equate=False),
        _Option(["-EvN", "EvN"],
                    "number of sequences for Viterbi Gumbel mu fit  [200]  (n>0)",
                    checker_function=lambda x: isinstance(x, int) and x>0,
                    equate=False),
        _Option(["-EfL", "EfL"],
                    "length of sequences for Forward exp tail tau fit  [100]  (n>0)",
                    checker_function=lambda x: isinstance(x, int) and x>0,
                    equate=False),
        _Option(["-EfN", "EfN"],
                    "number of sequences for Forward exp tail tau fit  [200]  (n>0)",
                    checker_function=lambda x: isinstance(x, int) and x>0,
                    equate=False),
        _Option(["-Eft", "Eft"],
                    "tail mass for Forward exponential tail tau fit  [0.04]  (0<x<1)",
                    checker_function=lambda x: isinstance(x, int) and 0<x<1,
                    equate=False),
        #Other options:
        _Switch(["--stall", "stall"],
                "arrest after start: for attaching debugger to process"),
        _Option(["-informat", "informat"],
                    "assert input alifile is in format <s> (no autodetect)",
                    equate=False),
        _Option(["-w_beta", "w_beta"],
                    "tail mass at which window length is determined",
                    equate=False),
        _Option(["-w_length", "w_length"],
                    "window length",
                    equate=False),
        _Option(["-maxinsertlen", "maxinsertlen"],
                    "pretend all inserts are length <= <n>",
                    equate=False),
        #Arguments
        _Argument(['out'],
                  "The HMM profile file to create",
                  filename=True, is_required=True),
        _Argument(['input'],
                  "An input multiple sequence alignment file",
                  filename=True, is_required=True)
        ]
        super(HMMBuildCommandline, self).__init__(cmd, **kwargs)
