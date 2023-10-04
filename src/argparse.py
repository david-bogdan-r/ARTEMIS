import argparse
import os
import multiprocessing as mp

def argParse(args):
    args = unify(args)
    args = parser.parse_args(args)
    setDefaultFormats(args)
    return args


def unify(args:'list'):
    ans = args.copy()
    i = 0
    while i < len(ans):
        v = ans[i]
        if '=' in v:
            k, v = v.split('=')
            if len(k) == 1:
                k = '-' + k
            else:
                k = '--' + k
            ans[i:i+1] = k, v
            i += 1

        i += 1

    return ans


def setDefaultFormats(args:'argparse.Namespace'):

    if args.rformat is None:
        fmt = os.path.splitext(args.r)[-1]
        if fmt in {'.cif', '.pdb', '.mmcif'}:
            args.rformat = fmt
        else:
            args.rformat = '.pdb'

    if args.qformat is None:
        fmt = os.path.splitext(args.q)[-1]
        if fmt in {'.cif', '.pdb', '.mmcif'}:
            args.qformat = fmt
        else:
            args.qformat = '.pdb'

    if args.saveformat is None:
        args.saveformat = args.qformat


    if args.saveres is None:
        if args.qresneg is None:
            args.saveres = args.qres
        else:
            args.saveres = args.qresneg

    if args.verbose and not args.permutation:
        args.permutation = True


def isfilepath(path:'str') -> 'str':
    if os.path.isfile(path):
        return path
    else:
        raise FileNotFoundError(path)

def fileformat(fmt:'str'):
    FMT = fmt.strip('.').upper()

    if FMT in {'PDB', 'CIF', 'MMCIF'}:
        if FMT == 'PDB':
            return '.pdb'
        else:
            return '.cif'
    else:
        raise argparse.ArgumentTypeError(
            "invalid choice: '{}' (choose from 'PDB', 'CIF', 'MMCIF')"
            .format(fmt)
        )

def threads(val):
    val = int(val)

    if val >= mp.cpu_count():
        return mp.cpu_count()
    elif val < 1:
        return mp.cpu_count()
    else:
        return val


parser = argparse.ArgumentParser(
    prog='ARTEMIS',
    description='''
    ARTEMIS (Using ARTEM to Infer Sequence alignment) 
    - is a tool for the RNA 3D structure superposition
    and determination of the structure-based sequence alignment.
    ''',
    add_help=False,
    usage= '''
    python3 artemis.py -r FILENAME -q FILENAME [OPTIONS]
    python3 artemis.py r=FILENAME q=FILENAME [OPTIONS]
    ''',
    epilog= '*'*os.get_terminal_size()[0] + '\n' + '''

    ARTEMIS uses a ChimeraX-like format to specify the residues of interest 
    using the "res" parameters:

    [#[INT]][/[STRING]][:[STRING][_INT[CHAR|_INT]] - The structure specification
                                                     format. The "res" 
                                                     parameters can be defined 
                                                     with a number 
                                                     of specifications 
                                                     separated by spaces and 
                                                     enclosed in double quotes.

        #[INT]                    == Model number
        /[STRING]                 == Chain identifier
        :[STRING][_INT[CHAR|_INT] == Residue(s) specification:
            
            :STRING     == Residue type    
            :_INT[CHAR] == Residue number [with insertion code]
            :_INT_INT   == Range of residue numbers
            
    Structure specification examples:

    rres="#1/B:_-10_20 #1/A"    - Consider the entire chain A from model 1 
                                  and the range of chain B residues with 
                                  numbers from -10 to 20 from model 1 as 
                                  the reference structure.
    qres="#"                    - Consider all the residues from all 
                                  the models in the "q" file as 
                                  the query structure.
    saveres="/C:_10_20 /C:_20A" - Save the chain C residues with numbers 
                                  from 10 to 20 and the chain C residue 
                                  with number 20A (A is the insertion code).
    rseed=:A                    - Use only the model 1 adenosines as the 
                                  single-residue seeds from the reference 
                                  structure. 


    Source: https://github.com/david-bogdan-r/ARTEMIS.git
    ''',
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument(
    '-h', '--h', '-H', '--H', '-help', '--help', 
    action='help',
    default=argparse.SUPPRESS,
    help='''
    Show this help message and exit.
    ''',
)

parser.add_argument(
    '-r', '--reference',
    metavar='REFERENCE',
    dest='r',
    required=True,
    type=isfilepath,
    help='''
    Path to a reference structure in PDB/mmCIF format. For faster 
    performance, it\'s advised to specify the largest of the two 
    structures as the reference structure.
    '''
)

parser.add_argument(
    '-q', '--query',
    metavar='QUERY',
    dest='q',
    required=True,
    type=isfilepath,
    help='''
    Path to a query structure, the one that ARTEMIS superimposes to
    the reference, in PDB/mmCIF format.
    '''
)

parser.add_argument(
    '-v', '--verbose',
    default=False,
    action='store_true',
    dest='verbose',
    required=False,
    help='''
    (default: verbose=False)
    Verbose mode.
    '''
)

parser.add_argument(
    '-p', '--permutation',
    default=False,
    action='store_true',
    dest='permutation',
    required=False,
    help='''
    (default: permutation=False)
    Permutation mode.
    If specified, ARTEMIS will add the permutation alignment details
    to the standard output, and save its output files
    to the saveto folder (if the folder is specified). 
    The mode is automatically activated if the TM-score of the query 
    structure for alignment with permutations is at least 10 percent higher 
    than the TM-score of the permutation-free alignment.
    '''
)

parser.add_argument(
    '-t', '--threads',
    default=mp.cpu_count(),
    dest='threads',
    required=False,
    type=threads,
    help='''
    (default: threads=%(default)d [CPU count])
    Number of CPUs to use.
    '''
)

parser.add_argument(
    '-m', '--matchrange',
    default=3.5,
    dest='matchrange',
    required=False,
    type=float,
    help='''
    (default: matchrange=3.5)
    The threshold used for searching the mutually closest residues.
    Only those pairs of residues that have their C3' atoms at a distance
    under the specified matchrange value can be added to the subset
    of the mutually closest residues. The higher matchrange value
    will produce more "noisy" matchings but won\'t miss anything. The lower
    matchrange value will produce more "clean" matchings but
    can miss something.
    ''',
)

parser.add_argument(
   '-d', '--stepdiv',
    default=None,
    dest='stepdiv',
    required=False,
    type=float,
    help='''
    (default: stepdiv=0 if QUERY length is less than 500nt else stepdiv=100)
    The step divider parameter. The parameter is used to speed up 
    the procedure for large structures. If stepdiv > 0, ARTEMIS
    will consider only each Sth reference residue as matching seed,
    where S = 1 + len(qres)//stepdiv. If the size of the query 
    structure exceeds 500 residues, stepdiv will be set to 100 
    by default.
    '''
)

parser.add_argument(
    '-n', '--nlargest',
    default=None,
    dest='nlargest',
    required=False,
    type=int,
    help='''
    (default: nlargest=(length of QUERY) if length of QUERY
    less than 500nt else nlargest=2*THREADS)
    Number of largest mutually closest residue sets for which 
    alignments are constructed.
    '''
)

parser.add_argument(
    '-s', '--shift',
    default=None,
    dest='shift',
    required=False,
    type=float,
    help='''
    (default: shift=3 if length of QUERY less than 500nt else shift=20)
    The shift value for the ScoreMatrix used in the Needleman-Wunsch
    algorithm. Larger shift provides higher coverage.
    '''
)

parser.add_argument(
    '-rformat', '--rformat',
    type=fileformat,
    help='''
    (default: rformat=[r extension]|PDB) 
    See -qformat.
    '''
)

parser.add_argument(
    '-qformat', '--qformat',
    default=None,
    type=fileformat,
    help='''
    (default: qformat=[q extension]|PDB) 
    The specification of the input coordinate file formats
    (case-insensitive). By default, ARTEMIS tries to infer the format
    from the extensions of the input filenames. ".pdb", ".cif",
    and ".mmcif" formats can be recognized (case-insensitive). In the case
    of any other extension ARTEMIS will treat the file as the PDB-format
    file by default. If the "rformat" ("qformat") parameter is specified
    and it\'s either "PDB", "CIF", or "MMCIF" (case-insensitive),
    ARTEMIS will treat the reference (query) coordinate file
    as the specified format.
    '''
)

parser.add_argument(
    '-rres', '--rres',
    default='#1',
    type=str,
    help='''
    (default: rres=#1)
    See -qres.
    '''
)

parser.add_argument(
    '-qres', '--qres',
    default='#1',
    type=str,
    help='''
    (default: qres=#1)
    The specification of the input reference (rres) and query (qres)
    structures. Only the specified residues will be considered as part
    of the structure and all the other residues will be ignored.
    See the format description at the end of the OPTIONS section.
    '''
)

parser.add_argument(
    '-rresneg', '--rresneg',
    default=None,
    type=str,
    help='''
    (default: rresneg=None)
    See -qresneg.
    '''
)

parser.add_argument(
    '-qresneg', '--qresneg',
    default=None,
    type=str,
    help='''
    (default: qresneg=None)
    The specification of the input reference (rresneg) and query (qresneg)
    structures. The specified residues will be ignored and all the other
    residues considered as part of the structure. If both "rres"
    and "rresneg" (or "qres" and "qresneg") are specified simultaneusly,
    ARTEMIS will ignore "rres" ("qres") and consider only "rresneg" ("qresneg").
    See the format description at the end ofthe OPTIONS section.
    '''
)

parser.add_argument(
    '-rseed', '--rseed',
    default=None,
    type=str,
    help='''
    (default: rseed=RRES|RRESNEG)
    See -qseed.
    '''
)

parser.add_argument(
    '-qseed', '--qseed',
    default=None,
    type=str,
    help='''
    (default: qseed=QRES|QRESNEG)
    The specification of the reference and query residues that ARTEMIS
    can use for single-residue matching seeds.
    See the format description at the end of the OPTIONS section.
    '''
)


parser.add_argument(
    '-saveto', '--saveto',
    default=None,
    type=str,
    help='''
    (default: saveto=None)
    Path to the output folder to save the coordinate files
    of the superimposed query structures along with the mutually
    closest residue subsets. If the specified folder does not exist,
    ARTEMIS will create it. If the folder is not specified,
    nothing will be saved.
    '''
)

parser.add_argument(
    '-saveformat', '--saveformat',
    default=None,
    type=fileformat,
    help='''
    (default: saveformat=QFORMAT)
    The specification of the format of the output coordinate files.
    By default, ARTEMIS will save the coordinate files in the same format
    as the query input file. If the "saveformat" parameter is specified
    and it\'s either "PDB", "CIF", or "MMCIF" (case-insensitive), ARTEMIS
    will save the output coordinate files in the specified format.
    '''
)

parser.add_argument(
    '-saveres', '--saveres',
    default=None,
    type=str,
    help='''
    (default: saveres=QRES|QRESNEG)
    The specification of the query structure residues that will be saved
    in the output coordinate files.
    See the format description at the end of the OPTIONS section.
    '''
)
