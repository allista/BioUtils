# coding=utf-8

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


class BlastDBcmdCommandline(AbstractCommandline):
    """Commandline object for the NCBI Blast+ blastdbcmd program.
    """
    def __init__(self, cmd="blastdbcmd", **kwargs):
        isstr = lambda x: isinstance(x, str)
        isint = lambda x: isinstance(x, int)
        self.parameters = [
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["-version", "version"],
                    "Print version number;  ignore other arguments."),
            # BLAST database options
            _Option(["-db", "db"],
                    "BLAST database name. Default = `nr'. Incompatible with: "
                    "list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path",
                    checker_function=isstr,
                    equate=False),
            _Option(["-dbtype", "dbtype"],
                    "<String, `guess', `nucl', `prot'>. Molecule type stored in BLAST database. Default = `guess'",
                    checker_function=isstr,
                    equate=False),
            # Retrieval options
            _Option(["-entry", "entry"],
                    "<String> Comma-delimited search string(s) of sequence identifiers: e.g.: 555, AC147927, 'gnl|dbname|tag', or 'all' to select all sequences in the database * Incompatible with:  entry_batch, pig, info, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path",
                    checker_function=isstr,
                    equate=False),
            _Option(["-entry_batch", "entry_batch"],
                    "<File_In> Input file for batch processing (Format: one entry per line, seq id  followed by optional space-delimited specifier(s) [range|strand|mask_algo_id] * Incompatible with:  entry, range, strand, mask_sequence_with, pig, info, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path",
                    filename=True,
                    equate=False),
            _Option(["-pig", "pig"],
                    "<Integer, >=0> PIG to retrieve * Incompatible with:  entry, entry_batch, target_only, info, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path",
                    checker_function=isint,
                    equate=False),
            _Switch(["-info", "info"],
                    "Print BLAST database information * Incompatible with:  entry, entry_batch, outfmt, strand, target_only, ctrl_a, get_dups, pig, range, mask_sequence, list, remove_redundant_dbs, recursive, list_outfmt, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path"),
            # Sequence retrieval configuration options
            _Option(["-range", "range"],
                    "<String> Range of sequence to extract in 1-based offsets (Format: start-stop, for start to end of sequence use start - ) * Incompatible with:  entry_batch, info, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path",
                    checker_function=isstr,
                    equate=False),
            _Option(["-strand", "strand"],
                    "<String, `minus', `plus'> Strand of nucleotide sequence to extract Default = `plus' * Incompatible with:  entry_batch, info, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path",
                    checker_function=isstr,
                    equate=False),
            _Option(["-mask_sequence_with", "mask_sequence_with"],
                    "<Integer> Produce lower-case masked FASTA using the algorithm ID specified * Incompatible with:  entry_batch",
                    checker_function=isint,
                    equate=False),
            # Output configuration options
            _Option(["-out", "out"],
                    "Output file for alignment. Default = `-'",
                    filename=True,
                    equate=False),
            _Option(["-outfmt", "outfmt"],
                    """Output format, where the available format specifiers are:
    %f means sequence in FASTA format
    %s means sequence data (without defline)
    %a means accession
    %g means gi
    %o means ordinal id (OID)
    %i means sequence id
    %t means sequence title
    %l means sequence length
    %h means sequence hash value
    %T means taxid
    %e means membership integer
    %L means common taxonomic name
    %S means scientific name
    %P means PIG
    %m means sequence masking data.
       Masking data will be displayed as a series of 'N-M' values
       separated by ';' or the word 'none' if none are available.
If '%f' is specified, all other format specifiers are ignored.
For every format except '%f', each line of output will correspond
to a sequence.
Default = `%f'
* Incompatible with:  info, list, recursive, remove_redundant_dbs,
list_outfmt, show_blastdb_search_path""",
                    checker_function=isstr,
                    equate=False),
            _Switch(["-target_only", "target_only"],
                    "Definition line should contain target entry only * Incompatible with:  pig, info, get_dups, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path"),
            _Switch(["-get_dups", "get_dups"],
                    "Retrieve duplicate accessions * Incompatible with:  info, target_only, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path"),
            # Output configuration options for FASTA format
            _Option(["-line_length", "line_length"],
                    "<Integer, >=1> Line length for output Default = `80' * Incompatible with:  list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path",
                    checker_function=isint,
                    equate=False),
            _Switch(["-ctrl_a", "ctrl_a"],
                    "Use Ctrl-A as the non-redundant defline separator * Incompatible with:  info, list, recursive, remove_redundant_dbs, list_outfmt, show_blastdb_search_path"),
            # BLAST database configuration and discovery options
            _Switch(["-show_blastdb_search_path", "show_blastdb_search_path"],
                    "Displays the default BLAST database search paths * Incompatible with:  entry, entry_batch, outfmt, strand, target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence, line_length, list, recursive, list_outfmt, remove_redundant_dbs"),
            _Option(["-list", "list"],
                    "List BLAST databases in the specified directory * Incompatible with:  info, entry, entry_batch, outfmt, strand, target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence, line_length, show_blastdb_search_path",
                    checker_function=isstr,
                    equate=False),
            _Switch(["-remove_redundant_dbs", "remove_redundant_dbs"],
                    "Remove the databases that are referenced by another alias file in the directory in question * Incompatible with:  info, entry, entry_batch, outfmt, strand, target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence, line_length, show_blastdb_search_path"),
            _Switch(["-recursive", "recursive"],
                    "Recursively traverse the directory structure to list available BLAST databases * Incompatible with:  info, entry, entry_batch, outfmt, strand, target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence, line_length, show_blastdb_search_path"),
            _Option(["-list_outfmt", "list_outfmt"],
                    """<String> Output format for the list option, where the available format specifiers
are:
    %f means the BLAST database absolute file name path
    %p means the BLAST database molecule type
    %t means the BLAST database title
    %d means the date of last update of the BLAST database
    %l means the number of bases/residues in the BLAST database
    %n means the number of sequences in the BLAST database
    %U means the number of bytes used by the BLAST database
For every format each line of output will correspond to a BLAST database.
Default = `%f %p'
* Incompatible with:  info, entry, entry_batch, outfmt, strand,
target_only, ctrl_a, get_dups, pig, range, db, info, mask_sequence,
line_length, show_blastdb_search_path""",
                    checker_function=isstr,
                    equate=False),
            _Option(["-logfile", "logfile"],
                    "<File_Out> File to which the program log should be redirected",
                    filename=True,
                    equate=False),
        ]
        super(BlastDBcmdCommandline, self).__init__(cmd, **kwargs)
