import argparse
import os
import multiprocessing as mp


def argParse(args):
    args = unify(args)
    args = argparser.parse_args(args)
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


argparser = argparse.ArgumentParser(
    prog='ARTEMIS',
    description='ARTEMIS (Using ARTEM to Infer Sequence alignment) - '
                'is a tool for the RNA 3D structure '
                'superposition and determination of '
                'the structure-based sequence alignment.',
    add_help=False,
    usage= '\n'
    'python3 artemis.py -r FILENAME -q FILENAME [OPTIONS]\n'
    'python3 artemis.py r=FILENAME q=FILENAME [OPTIONS]\n'
)

argparser.add_argument(
    '-h', '--h', '-H', '--H', '-help', '--help', 
    action='help',
    default=argparse.SUPPRESS,
    help='Show this help message and exit.'
)

argparser.add_argument(
    '-r', '--reference',
    metavar='',
    dest='r',
    required=True,
    type=isfilepath,
    help='Path to a reference structure in PDB/mmCIF format. For faster '
         'performance, it\'s advised to specify the largest of the two '
         'structures as the reference structure.'
)

argparser.add_argument(
    '-q', '--query',
    metavar='',
    dest='q',
    required=True,
    type=isfilepath,
    help='Path to a query structure, the one that ARTEMIS superimposes to '
         'the reference, in PDB/mmCIF format.'
)

argparser.add_argument(
    '-v', '--verbose',
    default=False,
    action='store_true',
    dest='verbose',
    required=False,
    help='(default: verbose=False) TODO'
)

argparser.add_argument(
    '-p', '--permutation',
    default=False,
    action='store_true',
    dest='permutation',
    required=False,
    help='(default: permutation=False) TODO'
)

argparser.add_argument(
    '-t', '--threads',
    metavar='',
    default=mp.cpu_count(),
    dest='threads',
    required=False,
    type=threads,
    help='(default: threads=%(default)d [CPU count]) Number of CPUs to use.'
)

argparser.add_argument(
    '-m', '--matchrange',
    metavar='',
    default=3.5,
    dest='matchrange',
    required=False,
    type=float,
    help='(default: matchrange=3.5) The threshold used for searching the mutually closest residues. '
        "Only those pairs of residues that have their C3' atoms at a distance "
        'under the specified matchrange value can be added to the subset '
        'of the mutually closest residues. The higher matchrange value '
        'will produce more "noisy" matchings but won\'t miss anything. The lower '
        'matchrange value will produce more "clean" matchings but '
        'can miss something.'
)

argparser.add_argument(
   '-d', '--stepdiv',
    metavar='',
    default=None,
    dest='stepdiv',
    required=False,
    type=float,
    help='(default: stepdiv=0 if QUERY length less than 500nt else stepdiv=100) To speed up the procedure of pairwise superpositions of structures '
        'to find sets of mutually closest residues, ARTEMIS can skip rseed '
        'residuals in steps of '
        '1 + (number of qres residues) // stepdiv . '
        'If divider is 0, ARTEMIS does not skip rseed residuals. '
        'If the number of qres residues exceeds 500 and stepdiv '
        'is not specified, then divider automatically becomes 100.'
)

argparser.add_argument(
    '-n', '--nlargest',
    metavar='',
    default=None,
    dest='nlargest',
    required=False,
    type=int,
    help='(default: nlargest=(length of QUERY) if length of QUERY less than 500nt else nlargest=2*THREADS) Number of largest mutually nearest sets of residues for which '
        'alignments are constructed to define best.'
)

argparser.add_argument(
    '-s', '--shift',
    metavar='',
    default=None,
    dest='shift',
    required=False,
    type=float,
    help='(default: shift=3 if length of QUERY less than 500nt else shift=20)The value by which the Score Matrix is shifted for Needleman-Wunsch. '
         'Larger shift, greater coverage.'
)

file = argparser.add_argument_group('input file handling')

file.add_argument(
    '-rformat', '--rformat',
    metavar='',
    type=fileformat,
    help='(default: rformat=[r extension]|PDB) See -qformat.'
)

file.add_argument(
    '-qformat', '--qformat',
    metavar='',
    default=None,
    type=fileformat,
    help='(default: qformat=[q extension]|PDB) The specification of the input coordinate file formats '
         '(case-insensitive). By default, ARTEMIS tries to infer the format '
         'from the extensions of the input filenames. ".pdb", ".cif", '
         'and ".mmcif" formats can be recognized (case-insensitive). In the case '
         'of any other extension ARTEMIS will treat the file as the PDB-format '
         'file by default. If the "rformat" ("qformat") parameter is specified '
         'and it\'s either "PDB", "CIF", or "MMCIF" (case-insensitive), '
         'ARTEMIS will treat the reference (query) coordinate file'
         'as the specified format.'
)

file.add_argument(
    '-rres', '--rres',
    metavar='',
    default='#1',
    type=str,
    help='(default: rres=#1) See -qres.'
)

file.add_argument(
    '-qres', '--qres',
    metavar='',
    default='#1',
    type=str,
    help='(default: qres=#1) The specification of the input reference (rres) and query (qres) '
         'structures. Only the specified residues will considered as part '
         'of the structure and all the other residues will be ignored. '
         'See the format description at the end of the OPTIONS section.'
)

file.add_argument(
    '-rresneg', '--rresneg',
    metavar='',
    default=None,
    type=str,
    help='(default: rresneg=None) See -qresneg.'
)

file.add_argument(
    '-qresneg', '--qresneg',
    metavar='',
    default=None,
    type=str,
    help='(default: qresneg=None) The specification of the input reference (rresneg) and query (qresneg) '
         'structures. The specified residues will be ignored and all the other '
         'residues considered as part of the structure. If both "rres" '
         'and "rresneg" (or "qres" and "qresneg") are specified simultaneusly, '
         'ARTEMIS will ignore "rres" ("qres") and consider only "rresneg" '
         '("qresneg"). '
         'See the format description at the end ofthe OPTIONS section.'
)

file.add_argument(
    '-rseed', '--rseed',
    metavar='',
    default=None,
    type=str,
    help='(default: rseed=RRES|RRESNEG) See -qseed.'
)

file.add_argument(
    '-qseed', '--qseed',
    metavar='',
    default=None,
    type=str,
    help='(default: qseed=QRES|QRESNEG) The specification of the reference and query residues that ARTEMIS '
         'can use for single-residue matching seeds. '
         'See the format description at the end of the OPTIONS section.'
)


save = argparser.add_argument_group('save results')

save.add_argument(
    '-saveto', '--saveto',
    metavar='',
    default=None,
    type=str,
    help='(default: saveto=None) Path to the output folder to save the coordinate files '
        'of the superimposed query structures along with the mutually '
        'closest residue subsets. If the specified folder does not exist, '
        'ARTEMIS will create it. If the folder is not specified, '
        'nothing will be saved.'
)

save.add_argument(
    '-saveformat', '--saveformat',
    metavar='',
    default=None,
    type=fileformat,
    help='(default: saveformat=QFORMAT)The specification of the format of the output coordinate files. '\
        'By default, ARTEMIS will save the coordinate files in the same format '\
        'as the query input file. If the "saveformat" parameter is specified '\
        'and it\'s either "PDB", "CIF", or "MMCIF" (case-insensitive), ARTEMIS '\
        'will save the output coordinate files in the specified format.'
)

save.add_argument(
    '-saveres', '--saveres',
    metavar='',
    default=None,
    type=str,
    help='(default: saveres=QRES|QRESNEG)The specification of the query structure residues that will be saved '
         'in the output coordinate files. '
         'See the format description at the end of the OPTIONS section.'
)
