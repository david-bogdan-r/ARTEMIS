import itertools
import multiprocessing as mp
import os
import sys
from datetime import datetime as dt
from functools import partial
from heapq import nlargest as nlargestf
from typing import Iterable

from src.DataModel import DataModel, dataModel
from src.Hit import getHit, mutuallyClosest
from src.ResRepr import ResRepr
from src.Align import Impose, Align

BASEDIR  = os.path.dirname(__file__)
HELP     = {'--H', '-H', '--h', '-h', '--help', '-help'}
FORMAT   = {
    'PDB'  : '.pdb',
    'CIF'  : '.cif',
    'MMCIF': '.cif'
}


r:'str|None' = None
q:'str|None' = None

rformat:'str' = 'PDB'
qformat:'str' = 'PDB'

rres:'str' = '#1'
qres:'str' = '#1'

rresneg:'str|bool' = False
qresneg:'str|bool' = False

rseed:'str' = rres
qseed:'str' = qres

saveto:'str|None' = None
saveformat:'str'  = qformat
saveres:'str'     = '#'

threads:'int' = 1

matchrange:'float' = 3.5

resrepr:'str' = 'artemis'

step_divider:'int|float' = 0

nlargest:'int' = 0

shift:'float'  = 3

class ARTEMIS:

    def __init__(self, r:'DataModel', q:'DataModel') -> 'None':

        self.r = r
        self.q = q

        self.hit = partial(
            getHit, rt=r.t, qm=q.m, matchrange=matchrange
        )

        self.impose = partial(Impose, rL=r.L, qL=q.L, 
                              rd2=r.d0**2, qd2=q.d0**2)

        self.align   = partial(Align, shift=shift, 
                               rc=r.c, qc=q.c, rm=r.m, qm=q.m, 
                               rseq=r.seq, qseq=q.seq,
                               impose=self.impose)

        self.ans  = None
        self.ans1 = None
        self.ans2 = None

    def run(self):

        tic = dt.now()

        r = self.r
        q = self.q

        seed = self.get_seed(step_divider=step_divider)

        h = map(
            mutuallyClosest, 
            nlargestf(nlargest, self.get_hit(seed), key=len)
        )
        h = [hh for hh in h if hh.size >= 6]

        if threads > 1:
            ali = pool.starmap(self.align, h)

        else:
            align = self.align
            ali   = [align(hh[0], hh[1]) for hh in h]

        ans1 = max(ali, key=lambda x: sum(x[0][1])) # type: ignore
        ans2 = max(ali, key=lambda x: sum(x[1][1])) # type: ignore

        toc = dt.now()
        time_total = round((toc - tic).total_seconds(), 4)

        rAli, qAli = ans2[1][-1]
        Lali = sum(i != '-' and j != '-' for i, j in zip(rAli, qAli))


        self.ans1 = ans1[0]
        self.ans2 = ans2[1]
        self.ans  = {
            'PDB1'  : r.name,
            'PDB2'  : q.name,
            'L1'    : r.L,
            'L2'    : q.L,
            'Lali'  : Lali,
            'TM1'   : round(ans2[1][1][0], 4), # type: ignore
            'TM2'   : round(ans2[1][1][1], 4), # type: ignore
            'RMSD'  : round(ans2[1][2], 4), # type: ignore
            'pTM1'  : round(ans1[0][1][0], 4),
            'pTM2'  : round(ans1[0][1][1], 4),
            'pRMSD' : round(ans1[0][2], 4),
            'time_total': time_total,
            'Ali1'  : rAli,
            'Ali2'  : qAli,
        }


    def get_hit(self, seed:'Iterable'):

        hit = self.hit

        if threads > 1:
            hits = pool.starmap(hit, seed)

        else:
            hits = [hit(*s) for s in seed]

        return hits


    def get_seed(self, step_divider:'int|float'=0):

        r = self.r
        q = self.q

        if isinstance(step_divider, int):
            if step_divider == 0:
                seed = itertools.product(
                    r.seed.dropna(), 
                    q.seed.dropna()
                )

            else:
                step = 1 + q.L // step_divider
                seed = itertools.product(
                    r.seed[::step].dropna(), 
                    q.seed.dropna()
                )

        elif isinstance(step_divider, float):
            step = int(step_divider * q.L)
            seed = itertools.product(
                r.seed[::step].dropna(), 
                q.seed.dropna()
            )

        else:
            seed = itertools.product(
                r.seed.dropna(), 
                q.seed.dropna()
            )

        return seed



if __name__ == '__main__':

    # Argument parsing 

    if any(arg in HELP for arg in sys.argv[1:]):
        with open(BASEDIR + '/help.txt', 'r') as file:
            print(file.read())
            exit()


    kwargs = dict(arg.split('=') for arg in sys.argv[1:])

    r = kwargs.get('r', r)
    q = kwargs.get('q', q)

    # rformat
    if isinstance(r, str) and os.path.isfile(r):
        rname, ext = os.path.splitext(r)

        rname = rname.split(os.sep)[-1]
        ext   = ext.strip('.').upper()

        if 'rformat' not in kwargs:
            if ext in FORMAT:
                rformat = FORMAT[ext]

            else:
                rformat = FORMAT[rformat]

        else:
            rformat:'str' = kwargs['rformat']
            rformat = rformat.strip('.').upper()

            if rformat in FORMAT:
                rformat = FORMAT[rformat]
            
            else:
                raise ValueError('Unexpected value rformat={}'
                                 .format(kwargs['rformat']))

    else:
        raise ValueError('Expected path value for argument r')

    # qformat
    if isinstance(q, str) and os.path.isfile(q):
        qname, ext = os.path.splitext(q)

        qname = qname.split(os.sep)[-1]
        ext   = ext.strip('.').upper()

        if 'qformat' not in kwargs:
            if ext in FORMAT:
                qformat = FORMAT[ext]

            else:
                qformat = FORMAT[qformat]

        else:
            qformat:'str' = kwargs['qformat']
            qformat = qformat.strip('.').upper()

            if qformat in FORMAT:
                qformat = FORMAT[qformat]
            
            else:
                raise ValueError('Unexpected value rformat={}'
                                 .format(kwargs['qformat']))

    else:
        raise ValueError('Expected path value for argument q')

    rres = kwargs.get('rres', rres)
    qres = kwargs.get('qres', qres)

    rresneg = kwargs.get('rresneg', rresneg)
    if isinstance(rresneg, str):
        rres    = rresneg
        rresneg = True

    qresneg = kwargs.get('qresneg', qresneg)
    if isinstance(qresneg, str):
        qres    = qresneg
        qresneg = True

    rseed = kwargs.get('rseed', rres)
    qseed = kwargs.get('qseed', qres)

    saveto     = kwargs.get('saveto', saveto)
    saveformat = kwargs.get('saveformat', qformat)
    saveres    = kwargs.get('saveres', qres)

    threads = float(kwargs.get('threads', threads)) # type: ignore
    if threads < 1:
        threads = mp.cpu_count()
    else:
        threads = min(int(threads), mp.cpu_count())

    matchrange = float(kwargs.get('matchrange', matchrange))

    if 'nlargest' in kwargs:
        nlargest = int(kwargs['nlargest'])

    # Set multiprocess start method

    if threads > 1:
        if 'fork' in mp.get_all_start_methods():
            mp.set_start_method('fork')

        pool = mp.Pool(threads)


    # Creating a function for residue representation

    resRepr = ResRepr(resrepr)


    # Load models

    rmodel = dataModel(
        rname, rformat, r,
        rres, rseed, '#', rresneg,
        resRepr
    )

    qmodel = dataModel(
        qname, qformat, q,
        qres, qseed, saveres, qresneg,
        resRepr
    )


    if not nlargest:
        nlargest = qmodel.L

    # ARTEMIS
    artemis = ARTEMIS(rmodel, qmodel)
    artemis.run()

    print(artemis.ans)
