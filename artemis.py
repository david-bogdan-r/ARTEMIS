import itertools
import multiprocessing as mp
import os
import sys
from datetime import datetime as dt
from functools import partial
from heapq import nlargest as nlargestf
from typing import Iterable
import numpy as np

from src.DataModel import DataModel, dataModel, MCBI
from src.Hit import getHit, mutuallyClosest, hitFromAli
from src.ResRepr import ResRepr
from src.Align import Impose, Align

output = '''
 ********************************************************************
 * ARTEMIS (Version 20230828)                                       *
 * using ARTEM to Infer Sequence alignment                          *
 * Reference: TODO                                                  *
 * Please email comments and suggestions to dav.bog.rom@gmail.com   *
 ********************************************************************

Name of structure r: {r}
Name of structure q: {q} (to be superimposed onto structure r)
Length of structure r: {Lr} residues
Length of structure q: {Lq} residues

Aligned length= {Lali}, RMSD= {RMSD:6.2f}, Seq_ID=n_identical/n_aligned= {seqID:4.3f}
TM-score= {TMr:6.5f} (normalized by length of structure r: L={Lr}, d0={d0r:.2f})
TM-score= {TMq:6.5f} (normalized by length of structure q: L={Lq}, d0={d0q:.2f})
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
{rAli}
{dist}
{qAli}

Alignment with permutations:
Aligned length= {pLali}
TM-score= {pTMr:6.5f} (normalized by length of structure r)
TM-score= {pTMq:6.5f} (normalized by length of structure q)
RMSD= {pRMSD:6.2f}
Seq_ID=n_identical/n_aligned= {pseqID:4.3f}

#Total CPU time is {time_total:5.2f} seconds
'''

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

step_divider:'float'
nlargest:'int'
shift:'float'

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

        self.ans1:'list' = []
        self.ans2:'list' = []

    def run(self):

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

        self.ans1 = ans1[0]
        self.ans2 = ans2[1] # type: ignore
        self.time_total = time_total

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

        if step_divider:
            step = 1 + int(q.L // step_divider)
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


    def show(self):

        rmodel = self.r
        qmodel = self.q

        ans1 = self.ans1
        ans2 = self.ans2

        rAli, qAli = ans2[-1]
        Lali = sum(i != '-' and j != '-' for i, j in zip(rAli, qAli))

        h = hitFromAli(rAli, qAli)

        rm = rmodel.m[h[0]]
        a, b = ans2[0]
        qm = np.dot(qmodel.m[h[1]], a) + b

        d = np.sqrt(((rm - qm) ** 2).sum(axis=1))

        dist = []
        i = -1
        j = -1
        m = -1
        for k in range(len(rAli)):
            b1 = rAli[k] != '-'
            b2 = qAli[k] != '-'

            if b1:
                i += 1
            if b2:
                j += 1

            if b1 and b2:
                m += 1
                if d[m] < 5:
                    dist.append(':')
                else:
                    dist.append('.')
            else:
                dist.append(' ')
        dist = ''.join(dist)

        n_identical = sum(c1 == c2 for c1, c2 in zip(rAli, qAli))
        n_aligned   = Lali

        seqID = n_identical / n_aligned

        rpath = r
        qpath = q

        rchain = rmodel.i.get_level_values(1).unique()
        if len(rchain) == 1:
            rpath += ':' + rchain[0]

        qchain = qmodel.i.get_level_values(1).unique()
        if len(qchain) == 1:
            qpath += ':' + qchain[0]

        h = ans1[-1]
        p_n_aligned = len(h[0])

        ri = rmodel.i[h[0]].to_list()
        qi = qmodel.i[h[1]].to_list()

        p_n_identical = 0
        for rr, qr in zip(ri, qi):
            if rr[2] == qr[2]:
                p_n_identical += 1

        pseqID = p_n_identical / p_n_aligned

        toc = dt.now()

        ans  = {
            'r'     : rpath,
            'q'     : qpath,
            'Lr'    : rmodel.L,
            'Lq'    : qmodel.L,
            'Lali'  : Lali,
            'TMr'   : ans2[1][0],
            'TMq'   : ans2[1][1],
            'RMSD'  : ans2[2],
            'rAli'  : rAli,
            'dist'  : dist,
            'qAli'  : qAli,
            'seqID' : seqID,
            'pLali' : p_n_aligned,
            'pTMr'  : ans1[1][0],
            'pTMq'  : ans1[1][1],
            'pRMSD' : ans1[2],
            'pseqID': pseqID,
            'd0r'   : rmodel.d0,
            'd0q'   : qmodel.d0,
            'time_total': (toc - tic).total_seconds(),
        }

        print(output.format(**ans),end='')


    def save(self) -> 'None':

        r = self.r
        q = self.q

        ans1 = self.ans1
        ans2 = self.ans2

        if os.path.isdir(saveto): # type: ignore
            files = set(os.listdir(saveto))
        else:
            files = set()

        a, b = ans2[0]
        qq = q.copy()
        qq.coord = np.dot(qq.coord, a) + b
        qq.atom_site = (qq.atom_site
                        .set_index(MCBI).loc[q.saveres]
                        .reset_index()[q.atom_site.columns])

        fname = '{}_to_{}{}'.format(q.name, r.name, saveformat)
        i = 0
        while fname in files:
            fname = '{}_to_{}_({}){}'.format(q.name, r.name, i, saveformat)
            i += 1

        if saveformat == '.pdb':
            qq.to_pdb('{}/{}'.format(saveto, fname))
        else:
            qq.to_cif('{}/{}'.format(saveto, fname))

        h  = hitFromAli(*ans2[-1])

        rm = self.r.m[h[0]]
        qm = np.dot(q.m[h[1]], a) + b
        d = np.sqrt(((rm - qm) ** 2).sum(axis=1))

        ri = self.r.i[h[0]].to_list()
        qi = self.q.i[h[1]].to_list()

        fname = '{}_to_{}.tsv'.format(q.name, r.name)
        i = 0
        while fname in files:
            fname = '{}_to_{}_({}).tsv'.format(q.name, r.name, i)
            i += 1

        text = '{}\tdist\t{}\n'.format(self.r, self.q)
        for i, (ii, jj) in enumerate(zip(ri, qi)):
            text += "{}\t{:.2}\t{}\n".format('.'.join(list(map(str, ii))),
                                             d[i],
                                             '.'.join(list(map(str, jj))))

        with open('{}/{}'.format(saveto, fname), 'w') as file:
            file.write(text)


        a, b = ans1[0]
        qq = q.copy()
        qq.coord = np.dot(qq.coord, a) + b
        qq.atom_site = (qq.atom_site
                        .set_index(MCBI).loc[q.saveres]
                        .reset_index()[q.atom_site.columns])

        fname = '{}_to_{}_p{}'.format(q.name, r.name, saveformat)
        i = 0
        while fname in files:
            fname = '{}_to_{}_p_({}){}'.format(q.name, r.name, i, saveformat)
            i += 1

        if saveformat == '.pdb':
            qq.to_pdb('{}/{}'.format(saveto, fname))
        else:
            qq.to_cif('{}/{}'.format(saveto, fname))


        h  = ans1[-1]

        rm = self.r.m[h[0]]
        qm = np.dot(q.m[h[1]], a) + b
        d = np.sqrt(((rm - qm) ** 2).sum(axis=1))

        ri = self.r.i[h[0]].to_list()
        qi = self.q.i[h[1]].to_list()

        fname = '{}_to_{}_p.tsv'.format(q.name, r.name)
        i = 0
        while fname in files:
            fname = '{}_to_{}_p_({}).tsv'.format(q.name, r.name, i)
            i += 1

        text = '{}\tdist\t{}\n'.format(self.r, self.q)
        for i, (ii, jj) in enumerate(zip(ri, qi)):
            text += "{}\t{:.2}\t{}\n".format('.'.join(list(map(str, ii))),
                                             d[i],
                                             '.'.join(list(map(str, jj))))

        with open('{}/{}'.format(saveto, fname), 'w') as file:
            file.write(text)



if __name__ == '__main__':

    tic = dt.now()

    # Argument parsing 

    if any(arg in HELP for arg in sys.argv[1:]):
        with open(BASEDIR + '/help.txt', 'r') as file:
            print(file.read())
            exit()


    args   = sys.argv[1:]
    kwargs = {}
    i = 0
    while i != len(args):
        k = args[i]
        if k.startswith('-'):
            k = k.strip('-')
            i += 1
            v = args[i]
            kwargs[k] = v
        else:
            k, v = k.split('=')
            kwargs[k] = v
        i += 1

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

    saveformat_ = saveformat.strip('.').upper()
    if saveformat_ in FORMAT:
        saveformat = FORMAT[saveformat_]
    else:
        saveformat = '.pdb'

    threads = float(kwargs.get('threads', threads)) # type: ignore
    if threads < 1:
        threads = mp.cpu_count()
    else:
        threads = min(int(threads), mp.cpu_count())

    matchrange = float(kwargs.get('matchrange', matchrange))


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


    default_state = qmodel.seed.notna().sum() >= 500

    nlargest     = int(kwargs.get('nlargest', (qmodel.L, 2*threads)[default_state]))
    shift        = float(kwargs.get('shift', (3, 20)[default_state]))
    step_divider = float(kwargs.get('step_divider', (0, 100)[default_state]))

    if qmodel.L >= 500:
        shift = 20


    # ARTEMIS

    artemis = ARTEMIS(rmodel, qmodel)
    artemis.run()

    if saveto:
        saveto = saveto.strip(os.sep)
        artemis.save()

    artemis.show()
