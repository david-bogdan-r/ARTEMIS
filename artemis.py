import datetime
import itertools
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from src import PDB
from src.artmsNW import GapGlobalAlign, OutputGapAlign1, OutputGapAlign2
from src.supos import MCS, RMSD, transform

BASEDIR = Path(__file__).parent.as_posix()
BASEFMT = {'PDB', 'CIF', 'MMCIF'}
HELPARG = {'--H', '-H', '--h', '-h', '--help', '-help'}


class Artemis:

    def __init__(self, *args, **kwargs):

        if args:
            if any(a in HELPARG for a in args):
                self.help()
                exit()
            kwargs |= dict(a.split('=') for a in args)

        threads = int(kwargs.get('threads', 1))
        if threads < 1 or threads > mp.cpu_count():
            threads = mp.cpu_count()
        if threads > 1:
            mp.set_start_method('fork')
        kwargs['threads'] = threads

        kwargs['rmsdmin'] = float(kwargs.get('rmsdmin', 0.))
        kwargs['rmsdmax'] = float(kwargs.get('rmsdmax', 1e10))

        kwargs['sizemin'] = float(kwargs.get('sizemin', 1.))
        kwargs['sizemax'] = float(kwargs.get('sizemax', 1e10))

        kwargs['rmsdsizemin'] = float(kwargs.get('rmsdsizemin', 0.))
        kwargs['rmsdsizemax'] = float(kwargs.get('rmsdsizemax', 1e10))

        kwargs['trim'] = bool(kwargs.get('trim', True))

    def help(self, path:'str' = BASEDIR + '/help.txt'):
        with open(path, 'r') as file:
            print(file.read())

if __name__ == '__main__':
    artmis = Artemis(*sys.argv[1:])
    print('some text')

    # args = sys.argv[1:]

    # if any(a in HELPARG for a in args):
    #     help()
    #     exit()

    # kwargs = dict(a.split('=') for a in args)

    # threads = int(kwargs.get('threads', 1))
    # if threads < 1 or threads > mp.cpu_count():
    #     threads = mp.cpu_count()
    # if threads > 1:
    #     mp.set_start_method('fork')

    # rmsdmin = float(kwargs.get('rmsdmin', 0.))
    # rmsdmax = float(kwargs.get('rmsdmax', 1e10))

    # sizemin = float(kwargs.get('sizemin', 1.))
    # sizemax = float(kwargs.get('sizemax', 1e10))

    # rmsdsizemin = float(kwargs.get('rmsdsizemin', 0.))
    # rmsdsizemax = float(kwargs.get('rmsdsizemax', 1e10))

    # trim = bool(kwargs.get('trim', True))

    # matchrange = float(kwargs.get('matchrange', 3.))

    # threshold = float(kwargs.get('threshold', 0.5))

    # weight = float(kwargs.get('weight', 0.))

    # k = kwargs.get('repr', '5')
    # RESREPR = PDB.ResRepr(k)

    # staronly = kwargs.get('staronly', 1)
