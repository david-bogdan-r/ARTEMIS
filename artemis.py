import glob
import itertools
import multiprocessing as mp
import os
import sys
import traceback
from functools import partial
from heapq import nlargest
from time import time
from typing import Callable, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist

try:
    from src.Kabsch import transform
    from src.NW import globalAlign
    from src.PDBio import BaseModel, getResSpec
    from src.resrepr import load_resrepr, resrepr
except:
    from .src.Kabsch import transform
    from .src.NW import globalAlign
    from .src.PDBio import BaseModel, getResSpec
    from .src.resrepr import load_resrepr, resrepr


SEEDPOOL    = 100_000
PMATCHRANGE = 8

BASEDIR  = os.path.dirname(__file__)

RESREPRSOURCE = [
    os.path.join(BASEDIR, 'src/resrepr/artemis_1.json'), # 3-atom representation of residues
    os.path.join(BASEDIR, 'src/resrepr/artemis_2.json'), # C3'-atom representation of residues
]
RESREPR      = [partial(resrepr, **d)
                for d in load_resrepr(*RESREPRSOURCE)]
RESTYPES     = [*set.union(*map(lambda x: set(x.keywords.keys()), RESREPR))]

MCBI = [
    'pdbx_PDB_model_num',
    'auth_asym_id',
    'auth_comp_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code'
]

INDEX = '{head}{config}{alignment}{distance_1}{permutation}{distance_2}{time}'

try:
    TERMINAL_WIDTH = os.get_terminal_size()[0]
except:
    TERMINAL_WIDTH = 79

HEAD  = '''
 ********************************************************************
 * ARTEMIS (Version 1.51)                                           *
 * using ARTEM to Infer Sequence alignment                          *
 * Reference: 10.1093/nar/gkae758                                   *
 * Please email comments and suggestions to dav.bog.rom@gmail.com   *
 ********************************************************************
'''

CONFIG = '''
Configuration:

r={r}
rformat={rformat}
rres={rres}
rresneg={rresneg}
rseed={rseed}

q={q}
qformat={qformat}
qres={qres}
qresneg={qresneg}
qseed={qseed}

matchrange={matchrange}
toplargest={toplargest}
shift={shift}
stepdiv={stepdiv}

threads={threads}
''' + '_'*TERMINAL_WIDTH + '\n'

ALIGNMENT = '''
Name of structure r: {rName}:{rChain}
Name of structure q: {qName}:{qChain} (to be superimposed onto structure r)
Length of structure r: {rLength} residues
Length of structure q: {qLength} residues
{dssr}
Aligned length= {aliLength}, RMSD= {RMSD:6.2f}, Seq_ID=n_identical/n_aligned= {Seq_ID:4.3f}
TM-score= {rTMscore:6.5f} (normalized by length of structure r: L={rLength}, d0={r_d0:.2f})
TM-score= {qTMscore:6.5f} (normalized by length of structure q: L={qLength}, d0={q_d0:.2f})

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
{rAlignment}
{distances}
{qAlignment}
'''

PERMUTATION = '\n' + '_'*TERMINAL_WIDTH + '''

Alignment with permutations:

Aligned length= {p_aliLength}
TM-score= {p_rTMscore:6.5f} (normalized by length of structure r)
TM-score= {p_qTMscore:6.5f} (normalized by length of structure q)
RMSD= {p_RMSD:6.2f}
Seq_ID=n_identical/n_aligned= {p_Seq_ID:4.3f}
'''

TIME = '''
#Total CPU time is {total_time:5.2f} seconds
'''


def hit(
    rx:'np.ndarray', qx:'np.ndarray',
    rt:'KDTree', qm:'np.ndarray',
    matchrange:'float'
) -> 'np.ndarray':

    a, b = transform(rx, qx)
    qt = KDTree(np.dot(qm, a) + b)

    h = rt.sparse_distance_matrix(
        qt,
        max_distance=matchrange,
        p=2,
        output_type='dict'
    )

    return mutuallyClosestHit(h)

def hitFromAli(ali1:'str', ali2:'str') -> 'np.ndarray':

    i = -1
    j = -1
    keys = []

    for c1, c2 in zip(ali1, ali2):

        b1 = c1 != '-'
        b2 = c2 != '-'

        if b1:
            i += 1
        if b2:
            j += 1

        if b1 and b2:
            keys.append([i, j])

    return np.array(keys).T

def mutuallyClosestHit(hit:'dict') -> 'np.ndarray':

    hk  = sorted(hit, key=hit.get) # type: ignore
    iset = set()
    jset = set()

    h = []
    for i, j in hk:
        if (i not in iset) and (j not in jset):
            iset.add(i)
            jset.add(j)
            h.append([i, j])

    h = np.array(h)

    return h

def impose(
    rm:'np.ndarray', qm:'np.ndarray',
    rL:'int', qL:'int',
    rd2:'float', qd2:'float',
    notmopt:'bool'
) -> 'dict':

    def diff(a, b) -> 'np.ndarray':

        d:'np.ndarray' = rm - (np.dot(qm, a) + b)

        return (d * d).sum(axis=1)

    def score(d2):

        rTM = (1 / (1 + d2 / rd2)).sum() / rL
        qTM = (1 / (1 + d2 / qd2)).sum() / qL

        RMSD = np.sqrt(np.sum(d2) / N)

        return (rTM, qTM), RMSD

    def opt(i:'int', j:'int'):

        n = j - i

        if n < 3:
            return [(None, None), (-1, -1), -1]

        nn = i + n // 2

        t = transform(rm[i:j], qm[i:j])
        d = diff(*t)
        s = score(d)

        m = d < qd2
        c = m.sum()

        if c > 2:
            tt = transform(rm[m], qm[m])
            ss = score(diff(*tt))

            if sum(ss[0]) > sum(s[0]):
                t = tt
                s = ss

        ans = max(
            [[t, *s], opt(i, nn), opt(nn, j)], 
            key=lambda x: sum(x[1]) # type: ignore
        )

        return ans

    N = len(rm)
    if not notmopt:
        o = opt(0, N)
    else:
        t = transform(rm, qm)
        d = diff(*t)
        s = score(d)
        o = [t, *s]

    ans = {
        'transform'  : o[0],
        'rTM'        : o[1][0],
        'qTM'        : o[1][1],
        'RMSD'       : o[2]
    }

    return ans

def align(
    ri:'np.ndarray' , qi:'np.ndarray',
    r:'DataModel'   , q:'DataModel',
    shift:'float',
    impose:'Callable'
) -> 'tuple[dict, dict]':

    rind = r.i[ri]
    qind = q.i[qi]

    rrepr1 = r.resrepr1.loc[rind]
    qrepr1 = q.resrepr1.loc[qind]

    rcar = rrepr1.notna().values
    qcar = qrepr1.notna().values

    car = rcar & qcar # type: ignore

    rx = np.vstack(rrepr1[car].values) # type: ignore
    qx = np.vstack(qrepr1[car].values) # type: ignore

    a, b = transform(rx, qx)

    rm = r.m
    qm = q.m

    dist = cdist(rm, np.dot(qm, a) + b)

    ri, qi = np.where(dist < PMATCHRANGE)

    h = mutuallyClosestHit(
        dict(zip(zip(ri, qi), dist[ri, qi]))
    ).T

    if not h.size:
        a = np.eye(3)
        b = np.zeros(3)

        ans1 = {
            'transform': (a, b),
            'rTM': -1,
            'qTM': -1,
            'RMSD': -1,
            'rAli': r.seq + '-'*len(q.seq),
            'qAli': q.seq + '-'*len(r.seq),
        }
        ans2 = {
            'transform': (a, b),
            'rTM': -1,
            'qTM': -1,
            'RMSD': -1,
            'rAli': r.i[0:0],
            'qAli': q.i[0:0]
        }
        return ans1, ans2

    ri, qi = h[:, h[0].argsort()]

    ans2 = impose(rm[ri], qm[qi])
    ans2['rAli'] = r.i[ri]
    ans2['qAli'] = q.i[qi]

    # scoremat  = -dist
    # scoremat -= scoremat.min()
    # scoremat -= min(scoremat[ri, qi]) - shift

    _shift   = max(dist[ri, qi]) + shift
    scoremat = -(dist - _shift)

    rAli, qAli = globalAlign(r.seq, q.seq, scoremat)
    h = hitFromAli(rAli, qAli)

    if h.size < 6:
        a = np.eye(3)
        b = np.zeros(3)

        ans1 = {
            'transform': (a, b),
            'rTM': -1,
            'qTM': -1,
            'RMSD': -1,
            'rAli': rAli,
            'qAli': qAli,
        }

    else:
        ri, qi = h

        ans1 = impose(rm[ri], qm[qi])
        ans1['rAli'] = rAli
        ans1['qAli'] = qAli

    return ans1, ans2


class DataModel(BaseModel):

    def __init__(self, path:'str', fmt:'str',
                 res:'str'='#1', resneg:'str|None'=None, seed:'str|None'=None
                 ) -> 'None':

        super().__init__(path, fmt)

        atom_site = self.atom_site

        resrepr1 = RESREPR[0](atom_site)
        resrepr2 = RESREPR[1](atom_site)

        mask = resrepr2.notna()
        resrepr2[mask] = resrepr2[mask].apply(lambda x: x.mean(axis=0))

        self.resrepr1 = resrepr1
        self.resrepr2 = resrepr2
        self._atom_site = atom_site[atom_site['auth_comp_id'].isin(RESTYPES)]

        self.set_res(res, resneg, seed)


    def set_res(self, res:'str'='#1', resneg:'str|None'=None, 
                seed:'str|None'=None)->'None':

        self.res    = res
        self.resneg = resneg

        resrepr1  = self.resrepr1
        resrepr2  = self.resrepr2
        atom_site = self._atom_site

        resi = (pd.MultiIndex
                .from_tuples(getResSpec(atom_site, res),
                             names=MCBI))

        if resneg:
            resi = resi.difference(
                pd.MultiIndex
                .from_tuples(getResSpec(atom_site, resneg),
                             names=MCBI),
                sort=False
            )

        if resi.empty:
            raise ValueError(
                'Empty residue set path={} (res={}, resneg={})'
                .format(self.path, res, resneg)
            )

        self.resi  = resi

        if seed is None:
            seedi = resi
            self.seed = seed
        else:
            seedi = getResSpec(atom_site, seed)
            seedi = pd.MultiIndex.from_tuples(seedi, names=MCBI)
            self.seed  = seed

        self.seedi = seedi

        s = resrepr1.loc[seedi]
        s = s[s.notna()].values
        self.s = s

        m = resrepr2.loc[resi]
        m = m[m.notna()]
        self.i = m.index

        m = np.vstack(m.values) # type: ignore
        self.m = m

        t = KDTree(m)
        self.t = t


    def __repr__(self) -> 'str':
        return '<{}{} DataModel>'.format(self, self.fmt)


    def __len__(self) -> 'int':
        return len(self.resi)


    def get_seq(self, repr2Exist:'bool'=True) -> 'str':

        if repr2Exist:

            seq = ''.join(
                [b if len(b) == 1 else 'M' 
                for b in self.i.get_level_values(2)]
            )

        else:
            seq = ''.join(
                [b if len(b) == 1 else 'M' 
                for b in self.resi.get_level_values(2)]
            )

        return seq


    def get_d0(self) -> 'float':

        L = len(self)

        if 30 <= L:
            return 0.6 * (L - 0.5)**(0.5) - 2.5

        elif 24 <= L < 30:
            return 0.7

        elif 20 <= L < 24:
            return 0.6

        elif 16 <= L < 20:
            return 0.5

        elif 12 <= L < 16:
            return 0.4

        else:
            return 0.3

    L   = property(__len__)
    seq = property(get_seq)
    d0  = property(get_d0)


class ARTEMIS:

    def __init__(self,
                 r:'DataModel|dict', q:'DataModel|dict', 
                 matchrange:'float'=3.5, threads:'int'=mp.cpu_count(),
                 toplargest:'int|None'=None, shift:'float|None'=None, 
                 stepdiv:'int|None'=None,
                 superonly:'bool'=False,
                 notmopt:'bool'=False,
                 ) -> 'None':

        self.tic        = time()

        self.matchrange = matchrange
        self.threads    = threads
        self.notmopt    = notmopt

        if threads > 1:
            self.pool = mp.Pool(threads)
        else:
            self.pool = None


        if isinstance(r, DataModel):
            self.r = r

        elif isinstance(r, dict):
            self.r = DataModel(**r)

        else:
            raise TypeError('r={}'.format(r))

        if isinstance(q, DataModel):
            self.q = q

        elif isinstance(q, dict):
            self.q = DataModel(**q)

        else:
            raise TypeError('q={}'.format(q))

        M = min(self.r.L, self.q.L)
        if M >= 500:

            if toplargest is None:
                toplargest = 2 * threads

            if shift is None:
                shift = 20

            if stepdiv is None:
                stepdiv = 100

        else:
            if toplargest is None:
                toplargest = M

            if shift is None:
                shift = 3

            if stepdiv is None:
                stepdiv = 0

        self.toplargest:'int' = toplargest # type: ignore
        self.shift    = shift
        self.stepdiv  = stepdiv

        self.hit = partial(
            hit,
            rt=self.r.t, qm=self.q.m, matchrange=self.matchrange # type: ignore
        )

        self.impose = partial(
            impose, 
            rL=self.r.L, qL=self.q.L,
            rd2=self.r.d0 ** 2, qd2=self.q.d0 ** 2,
            notmopt=self.notmopt
        )

        self.align = partial(
            align,
            r=self.r, q=self.q,
            shift=self.shift,
            impose=self.impose
        )

        self.ans1 = {}
        self.ans2 = {}

        self.perm = False
        self.superonly = superonly


    def insertNaNGap(self, rAli, qAli):

        r = self.r
        q = self.q

        rmsk = r.resi.isin(r.i)
        qmsk = q.resi.isin(q.i)

        newri = np.where(rmsk)[0]
        newqi = np.where(qmsk)[0]

        rali = []
        i = 0
        for c in rAli:
            if c == '-':
                rali.append(-1)
            else:
                rali.append(newri[i])
                i += 1

        qali = []
        i = 0
        for c in qAli:
            if c == '-':
                qali.append(-1)
            else:
                qali.append(newqi[i])
                i += 1

        i = 0
        j = 0
        while j != len(rali):
            ii = rali[j]
            if ii == -1:
                j += 1
                continue
            while i < ii:
                rali.insert(j, i)
                qali.insert(j, -1)
                j += 1
                i += 1
            else:
                i = ii+1
            j += 1

        i = 0
        j = 0
        while j != len(qali):
            ii = qali[j]
            if ii == -1:
                j += 1
                continue
            while i < ii:
                qali.insert(j, i)
                rali.insert(j, -1)
                j += 1
                i += 1
            else:
                i = ii+1
            j += 1

        rAli = ''
        for i in rali:
            if i == -1:
                rAli += '-'
            else:
                c = r.resi[i][2]
                if len(c) > 1:
                    if c in {'DA', 'DG', 'DC', 'DT'}:
                        c = c[-1]
                    else:
                        c = 'M'
                rAli += c

        qAli = ''
        for i in qali:
            if i == -1:
                qAli += '-'
            else:
                c = q.resi[i][2]
                if len(c) > 1:
                    if c in {'DA', 'DG', 'DC', 'DT'}:
                        c = c[-1]
                    else:
                        c = 'M'
                qAli += c

        # broken residue to lowercase
        rAli_lc = ''
        i = -1
        for c in rAli:
            if c != '-':
                i += 1

                if not rmsk[i]:
                    rAli_lc += c.lower()
                else:
                    rAli_lc += c

            else:
                rAli_lc += c

        qAli_lc = ''
        i = -1
        for c in qAli:
            if c != '-':
                i += 1

                if not qmsk[i]:
                    qAli_lc += c.lower()
                else:
                    qAli_lc += c

            else:
                qAli_lc += c

        return rAli_lc, qAli_lc

    def get_seed(self):

        if self.stepdiv == 0:
            rseed = self.r.s
        else:
            step = 1 + int(self.q.L // self.stepdiv)
            rseed = self.r.s[::step]

        seed = itertools.product(
            rseed, self.q.s
        )

        count = len(rseed) * len(self.q.s)

        return seed, count

    def get_hit(self, seed:'Iterable'):

        hit = self.hit

        if self.pool is not None:
            hits = self.pool.starmap(hit, seed)

        else:
            hits = [hit(*s) for s in seed]

        return hits

    def run(self):

        if self.superonly:

            r = self.r
            q = self.q

            n = len(r.m)
            if n != len(q.m):
                raise ValueError("Sequence lengths are not equal, disable superonly mode.")

            if n < 3:
                raise ValueError("Sequence lengths < 3, superposition is not possible.")

            ri = np.arange(n, dtype=int)
            qi = np.arange(n, dtype=int)

            ans2 = self.impose(r.m, q.m)
            ans2['rAli'] = r.i[ri]
            ans2['qAli'] = q.i[qi]


            ans1 = ans2.copy()
            ans1['rAli'] = r.seq
            ans1['qAli'] = q.seq
            rAli, qAli = self.insertNaNGap(ans1['rAli'], ans1['qAli'])
            ans1['rAli'] = rAli
            ans1['qAli'] = qAli

            self.ans1 = ans1
            self.ans2 = ans2
            return

        seed, count = self.get_seed()

        i = 0
        h = []

        while i + SEEDPOOL < count:
            step = [next(seed) for _ in range(SEEDPOOL)]
            h = nlargest(
                self.toplargest,
                h + self.get_hit(step),
                key=len
            )
            i += SEEDPOOL

        h = nlargest(
            self.toplargest,
            h + self.get_hit(seed),
            key=len
        )

        h = [hh.T for hh in h if len(hh) >= 3]

        if not h:
            raise Exception(
            'No correspondence between sets of residues was found. '\
            'Try other argument values.'
            )

        if self.pool is not None:
            ali = self.pool.starmap(self.align, h)
        else:
            align = self.align
            ali = [align(hh[0], hh[1]) for hh in h]

        ans1 = max(
            ali,
            key=lambda x: x[0]['rTM'] + x[0]['qTM']
        )[0]
        rAli, qAli = self.insertNaNGap(ans1['rAli'], ans1['qAli'])
        ans1['rAli'] = rAli
        ans1['qAli'] = qAli

        ans2 = max(
            ali,
            key=lambda x: x[1]['rTM'] + x[1]['qTM']
        )[1]

        if (ans2['qTM'] - ans1['qTM']) / abs(ans1['qTM']) >= 0.1:
            self.perm = True

        self.ans1 = ans1
        self.ans2 = ans2

    def get_lifetime(self) -> 'float':
        toc = time()
        return toc - self.tic

    lifetime = property(get_lifetime)

    def get_config(self):

        r = self.r
        q = self.q

        if r.seed is None:
            if r.resneg:
                r.seed = 'rres-rresneg(default)'
            else:
                r.seed = 'rres(default)'

        if q.seed is None:
            if q.resneg:
                q.seed = 'qres-qresneg(default)'
            else:
                q.seed = 'qres(default)'

        config = {
            'r'         : r.path,
            'rformat'   : r.fmt,
            'rres'      : r.res if r.res else '',
            'rresneg'   : r.resneg if r.resneg else '',
            'rseed'     : r.seed,

            'q'         : q.path,
            'qformat'   : q.fmt,
            'qres'      : q.res if q.res else '',
            'qresneg'   : q.resneg if q.resneg else '',
            'qseed'     : q.seed,

            'matchrange': self.matchrange,
            'toplargest': self.toplargest,
            'shift'     : self.shift,
            'stepdiv'   : self.stepdiv,

            'threads'   : self.threads,
        }

        return config


    def get_alignment(self) -> 'dict':

        r = self.r
        q = self.q

        ans1  = self.ans1
        if not ans1:
            self.run()
            ans1 = self.ans1

        rChain = ', '.join(
            r.resi
            .get_level_values(1)
            .unique()
            .astype(str)
        )

        qChain = ', '.join(
            q.resi
            .get_level_values(1)
            .unique()
            .astype(str)
        )

        rAli = ans1['rAli']
        qAli = ans1['qAli']

        aliLength = sum(i != '-' and j != '-' for i, j in zip(rAli, qAli))

        if aliLength:
            n_identical = sum(c1 == c2 for c1, c2 in zip(rAli, qAli))
            n_aligned   = aliLength
            Seq_ID = n_identical / n_aligned


            table = self.get_distance_1()
            distances = []
            i = 0
            for k in range(len(rAli)):
                if rAli[k] != '-' and qAli[k] != '-':
                    d = table['dist'][i]
                    if d < 5:
                        distances.append(':')
                    else:
                        distances.append('.')

                    i += 1
                else:
                    distances.append(' ')

        else:
            Seq_ID = 0
            distances = [' '] * len(rAli)

        distances = ''.join(distances)


        alignment = {
            'rName'     : r.path,
            'qName'     : q.path,
            'rChain'    : rChain,
            'qChain'    : qChain,
            'rLength'   : r.L,
            'qLength'   : q.L,
            'aliLength' : aliLength,
            'RMSD'      : ans1['RMSD'],
            'Seq_ID'    : Seq_ID,
            'rTMscore'  : ans1['rTM'],
            'r_d0'      : r.d0,
            'qTMscore'  : ans1['qTM'],
            'q_d0'      : q.d0,
            'rAlignment': rAli,
            'distances' : distances,
            'qAlignment': qAli,
            'transform' : ans1['transform']
        }

        return alignment

    def get_permutation(self) -> 'dict':

        ans2 = self.ans2

        if not ans2:
            self.run()
            ans2 = self.ans2

        p_aliLength = len(ans2['rAli'])
        p_rTMscore  = ans2['rTM']
        p_qTMscore  = ans2['qTM']
        p_RMSD      = ans2['RMSD']

        p_n_identical = 0
        for rc, qc in zip(ans2['rAli'], ans2['qAli']):
            if rc[2] == qc[2]:
                p_n_identical += 1
        p_Seq_ID = p_n_identical / p_aliLength

        permutation = {
            'p_aliLength'   : p_aliLength,
            'p_rTMscore'    : p_rTMscore,
            'p_qTMscore'    : p_qTMscore,
            'p_RMSD'        : p_RMSD,
            'p_Seq_ID'      : p_Seq_ID,
            'transform'     : ans2['transform']
        }

        return permutation

    def get_distance_1(self) -> 'pd.DataFrame':

        r = self.r
        q = self.q

        ans1  = self.ans1
        if not ans1:
            self.run()
            ans1 = self.ans1

        rAli = ans1['rAli']
        qAli = ans1['qAli']

        aliLength = sum(i != '-' and j != '-' for i, j in zip(rAli, qAli))

        rname = r.name
        qname = q.name
        if rname == qname:
                rname += '_reference'
                qname += '_query'

        if aliLength:
            a, b   = ans1['transform']
            ri, qi = hitFromAli(rAli, qAli)

            ri = r.resi[ri]
            qi = q.resi[qi]

            rcoord = np.vstack(r.resrepr2[ri].values) # type: ignore
            qcoord = np.vstack(q.resrepr2[qi].values) # type: ignore
            qcoord = np.dot(qcoord, a) + b

            d  = np.sqrt(((rcoord - qcoord) ** 2).sum(axis=1))
            ri = ['.'.join([*map(str, label)]) for label in ri]
            qi = ['.'.join([*map(str, label)]) for label in qi]

            table = pd.DataFrame(
                {
                    rname: ri,
                    'dist': d,
                    qname: qi
                }
            )

        else:
            table = pd.DataFrame(
                {
                    rname: [],
                    'dist': [],
                    qname: []
                }
            )

        table[rname] = table[rname].str.rstrip('?')
        table[qname] = table[qname].str.rstrip('?')

        return table

    def get_distance_2(self) -> 'pd.DataFrame':

        r = self.r
        q = self.q

        ans2 = self.ans2
        rAli = ans2['rAli']
        qAli = ans2['qAli']

        rname = r.name
        qname = q.name
        if rname == qname:
                rname += '_reference'
                qname += '_query'

        if len(rAli) != 0:

            a, b = ans2['transform']

            rcoord = np.vstack(r.resrepr2[rAli].values) # type: ignore
            qcoord = np.vstack(q.resrepr2[qAli].values) # type: ignore
            qcoord = np.dot(qcoord, a) + b

            d  = np.sqrt(((rcoord - qcoord) ** 2).sum(axis=1))
            ri = ['.'.join([*map(str, label)]) for label in rAli]
            qi = ['.'.join([*map(str, label)]) for label in qAli]

            table = pd.DataFrame(
                {
                    rname: ri,
                    'dist': d,
                    qname: qi
                }
            )

        else:
            table = pd.DataFrame(
                {
                    rname: [],
                    'dist': [],
                    qname: []
                }
            )

        table[rname] = table[rname].str.rstrip('?')
        table[qname] = table[qname].str.rstrip('?')

        return table

    def get_time(self) -> 'dict':

        total_time = {
            'total_time': self.lifetime
        }

        return total_time

    def show(self, verbose:'bool'=False, permutation:'bool'=False, dssr:'bool'=False) -> 'str':

        if dssr:
            dssr_s  = 'rDSSR: ' + ', '.join(['.'.join(map(str, i)) for i in self.r.resi]) + '\n'
            dssr_s += 'qDSSR: ' + ', '.join(['.'.join(map(str, i)) for i in self.q.resi]) + '\n\n'
        else:
            dssr_s = ''

        index = {
            'head'      : HEAD,
            'alignment' : ALIGNMENT.format(dssr=dssr_s, **self.get_alignment()),
            'time'      : TIME.format(total_time=self.lifetime),
        }

        if verbose:
            permutation=True

        if permutation:
            index['permutation'] = PERMUTATION.format(**self.get_permutation())
            index['distance_2'] = '\nDistance table:\n\n'\
                + (self
                    .get_distance_2()
                    .to_string(
                        index=False,
                        justify='center',
                        float_format='{:.2f}'.format
                        ) + '\n'
                    )
        else:
            index['permutation'] = ''
            index['distance_2']  = ''

        if verbose:
            index['config']     = CONFIG.format(**self.get_config())

            index['distance_1'] = '\nDistance table:\n\n'\
                + (self
                   .get_distance_1()
                   .to_string(
                       index=False,
                       justify='center',
                       float_format='{:.2f}'.format
                       )
                   )

        else:
            index['config']     = ''
            index['distance_1'] = ''

        return INDEX.format(**index)

    def save(self, 
             saveto     :'str'='.',
             saveres    :'str'='',  # type: ignore
             saveformat :'str'='',
             permutation:'bool'=False,
             ID=None) -> 'None':

        r = self.r
        q = self.q

        if not saveres:
            # if isinstance(q.resneg, str):
            #     saveres:'str' = q.resneg
            # else:
            #     saveres:'str' = q.res

            # saveresi = q.resi

            saveresi = pd.MultiIndex.from_tuples(
                getResSpec(q.atom_site, q.res)
            )
            if q.resneg:
                saveresi = saveresi.difference(
                    pd.MultiIndex.from_tuples(
                        getResSpec(q.atom_site, q.resneg)
                    )
                )

        else:
            saveresi = getResSpec(q.atom_site, saveres)

        if not saveformat:
            saveformat = q.fmt

        # saveresi = getResSpec(q.atom_site, saveres)

        ans1 = self.ans1
        ans2 = self.ans2

        rAli = ans1['rAli']
        qAli = ans1['qAli']

        h  = hitFromAli(rAli, qAli)

        if h.size:
            a, b = ans1['transform']

            qq       = q.copy()
            qq.coord = np.dot(qq.coord, a) + b
            qq.atom_site = (qq.atom_site
                            .set_index(MCBI).loc[saveresi]
                            .reset_index()[q.atom_site.columns])

            if ID is None:
                fname = '{}_to_{}{}'.format(q.name, r.name, saveformat)
            else:
                fname = '{}_to_{}_{}{}'.format(q.name, r.name, ID, saveformat)

            if saveformat == '.cif':
                qq.to_cif('{}/{}'.format(saveto, fname))
            else:
                qq.to_pdb('{}/{}'.format(saveto, fname))

            table = self.get_distance_1()

            if ID is None:
                fname = '{}_to_{}.tsv'.format(q.name, r.name)
            else:
                fname = '{}_to_{}_{}.tsv'.format(q.name, r.name, ID)

            table.to_csv(saveto + '/' + fname, sep='\t', 
                         float_format='{:.3f}'.format, index=False)

            if ID is None:
                fname =  '{}_to_{}.png'.format(q.name, r.name)
            else:
                fname =  '{}_to_{}_{}.png'.format(q.name, r.name, ID)

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot()
            ax.scatter(h[1] + 0.5, h[0] + 0.5,
                    marker='s',
                    s = max(518400 / (r.L * q.L), 0.2),
                    edgecolor='none',
                    c='black')
            ax.set_xlim(-1, q.L + 1)
            ax.set_ylim(-1, r.L + 1)
            fig.gca().invert_yaxis()
            ax.tick_params(
                top=True, 
                labeltop=True, 
                bottom=False, 
                labelbottom=False
            )
            ax.set_title(q.name, fontsize=18)
            ax.set_ylabel(r.name, fontsize=18)
            fig.savefig(saveto + '/' + fname, dpi=200, bbox_inches='tight')
            plt.close(fig)

        if permutation:
            a, b = ans2['transform']

            qq = q.copy()
            qq.coord = np.dot(qq.coord, a) + b
            qq.atom_site = (qq.atom_site
                            .set_index(MCBI).loc[saveresi]
                            .reset_index()[q.atom_site.columns])

            if ID is None:
                fname = '{}_to_{}_ti{}'.format(q.name, r.name, saveformat)
            else:
                fname = '{}_to_{}_{}_ti{}'.format(q.name, r.name, ID, saveformat)

            if saveformat == '.pdb':
                qq.to_pdb('{}/{}'.format(saveto, fname))
            else:
                qq.to_cif('{}/{}'.format(saveto, fname))

            table = self.get_distance_2()

            if ID is None:
                fname = '{}_to_{}_ti.tsv'.format(q.name, r.name)
            else:
                fname = '{}_to_{}_{}_ti.tsv'.format(q.name, r.name, ID)


            table.to_csv(saveto + '/' + fname, sep='\t',
                         float_format='{:.3f}'.format, index=False)

            if ID is None:
                fname =  '{}_to_{}_ti.png'.format(q.name, r.name)
            else:
                fname =  '{}_to_{}_{}_ti.png'.format(q.name, r.name, ID)

            fig = plt.figure(figsize=(10, 10))
            rAli, qAli = self.ans2['rAli'], self.ans2['qAli']
            lAli = len(rAli)
            rd = dict(zip(r.i, range(len(r.i))))
            qd = dict(zip(q.i, range(len(q.i))))
            ri, qi = [], []
            for k in range(lAli):
                ri.append(rd[rAli[k]])
                qi.append(qd[qAli[k]])
            ri = np.array(ri)
            qi = np.array(qi)
            # amat[ri, qi] = 1
            ax = fig.add_subplot()
            ax.scatter(qi+0.5, ri+0.5,
                       marker='s',
                       s = max(518400 / (r.L * q.L), 0.2),
                       edgecolor='none',
                       c='black')
            ax.set_xlim(-1, q.L + 1)
            ax.set_ylim(-1, r.L + 1)
            # ax.pcolor(amat, cmap='gray_r')
            fig.gca().invert_yaxis()
            ax.tick_params(
                top=True, 
                labeltop=True, 
                bottom=False, 
                labelbottom=False
            )
            ax.set_title(q.name, fontsize=18)
            ax.set_ylabel(r.name, fontsize=18)
            fig.savefig(saveto + '/' + fname, dpi=200, bbox_inches='tight')
            plt.close(fig)

    def save_seqord(self, 
                    saveto     :'str'='.',
                    saveres    :'str'='',  # type: ignore
                    saveformat :'str'='',
                    ID=None) -> 'None':
        r = self.r
        q = self.q

        if not saveres:
            # saveresi = q.resi
            saveresi = pd.MultiIndex.from_tuples(
                getResSpec(q.atom_site, q.res)
            )
            if q.resneg:
                saveresi = saveresi.difference(
                    pd.MultiIndex.from_tuples(
                        getResSpec(q.atom_site, q.resneg)
                    )
                )

        else:
            saveresi = getResSpec(q.atom_site, saveres)

        if not saveformat:
            saveformat = q.fmt

        ans1 = self.ans1

        rAli = ans1['rAli']
        qAli = ans1['qAli']

        h  = hitFromAli(rAli, qAli)

        if h.size:
            a, b = ans1['transform']

            qq       = q.copy()
            qq.coord = np.dot(qq.coord, a) + b
            qq.atom_site = (qq.atom_site
                            .set_index(MCBI).loc[saveresi]
                            .reset_index()[q.atom_site.columns])

            if ID is None:
                fname = '{}_to_{}{}'.format(q.name, r.name, saveformat)
            else:
                fname = '{}_to_{}_{}{}'.format(q.name, r.name, ID, saveformat)

            if saveformat == '.cif':
                qq.to_cif('{}/{}'.format(saveto, fname))
            else:
                qq.to_pdb('{}/{}'.format(saveto, fname))

            table = self.get_distance_1()

            if ID is None:
                fname = '{}_to_{}.tsv'.format(q.name, r.name)
            else:
                fname = '{}_to_{}_{}.tsv'.format(q.name, r.name, ID)

            table.to_csv(saveto + '/' + fname, sep='\t', 
                         float_format='{:.3f}'.format, index=False)

            if ID is None:
                fname =  '{}_to_{}.png'.format(q.name, r.name)
            else:
                fname =  '{}_to_{}_{}.png'.format(q.name, r.name, ID)

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot()
            ax.scatter(h[1] + 0.5, h[0] + 0.5,
                        marker='s',
                        edgecolor='none',
                        s = max(518400 / (r.L * q.L), 0.2),
                        c='black')
            ax.set_xlim(-1, q.L + 1)
            ax.set_ylim(-1, r.L + 1)
            fig.gca().invert_yaxis()
            ax.tick_params(
                top=True, 
                labeltop=True, 
                bottom=False, 
                labelbottom=False
            )
            ax.set_title(q.name, fontsize=18)
            ax.set_ylabel(r.name, fontsize=18)
            fig.savefig(saveto + '/' + fname, dpi=200, bbox_inches='tight')
            plt.close(fig)

    def save_topind(self, 
                    saveto     :'str'='.',
                    saveres    :'str'='',  # type: ignore
                    saveformat :'str'='',
                    ID=None) -> 'None':

        r = self.r
        q = self.q

        if not saveres:
            # saveresi = q.resi
            saveresi = pd.MultiIndex.from_tuples(
                getResSpec(q.atom_site, q.res)
            )
            if q.resneg:
                saveresi = saveresi.difference(
                    pd.MultiIndex.from_tuples(
                        getResSpec(q.atom_site, q.resneg)
                    )
                )
        else:
            saveresi = getResSpec(q.atom_site, saveres)

        if not saveformat:
            saveformat = q.fmt

        ans2 = self.ans2
        a, b = ans2['transform']
        qq = q.copy()
        qq.coord = np.dot(qq.coord, a) + b
        qq.atom_site = (qq.atom_site
                        .set_index(MCBI).loc[saveresi]
                        .reset_index()[q.atom_site.columns])

        if ID is None:
            fname = '{}_to_{}_ti{}'.format(q.name, r.name, saveformat)
        else:
            fname = '{}_to_{}_{}_ti{}'.format(q.name, r.name, ID, saveformat)

        if saveformat == '.pdb':
            qq.to_pdb('{}/{}'.format(saveto, fname))
        else:
            qq.to_cif('{}/{}'.format(saveto, fname))

        table = self.get_distance_2()

        if ID is None:
            fname = '{}_to_{}_ti.tsv'.format(q.name, r.name)
        else:
            fname = '{}_to_{}_{}_ti.tsv'.format(q.name, r.name, ID)


        table.to_csv(saveto + '/' + fname, sep='\t',
                        float_format='{:.3f}'.format, index=False)

        if ID is None:
            fname =  '{}_to_{}_ti.png'.format(q.name, r.name)
        else:
            fname =  '{}_to_{}_{}_ti.png'.format(q.name, r.name, ID)

        fig = plt.figure(figsize=(10, 10))
        rAli, qAli = self.ans2['rAli'], self.ans2['qAli']
        lAli = len(rAli)
        rd = dict(zip(r.i, range(len(r.i))))
        qd = dict(zip(q.i, range(len(q.i))))
        ri, qi = [], []
        for k in range(lAli):
            ri.append(rd[rAli[k]])
            qi.append(qd[qAli[k]])
        # amat[ri, qi] = 1
        ri = np.array(ri)
        qi = np.array(qi)
        ax = fig.add_subplot()
        # ax.pcolor(amat, cmap='gray_r')
        ax.scatter(
            qi+0.5, ri+0.5,
            marker='s',
            s = max(518400 / (r.L * q.L), 0.2),
            edgecolor='none',
            c='black'
        )
        ax.set_xlim(-1, q.L + 1)
        ax.set_ylim(-1, r.L + 1)
        fig.gca().invert_yaxis()
        ax.tick_params(
            top=True, 
            labeltop=True, 
            bottom=False, 
            labelbottom=False
        )
        ax.set_title(q.name, fontsize=18)
        ax.set_ylabel(r.name, fontsize=18)
        fig.savefig(saveto + '/' + fname, dpi=200, bbox_inches='tight')
        plt.close(fig)


def pathparser(mask:'str') -> 'list':

    paths = []

    bashlist = glob.glob(mask)
    if bashlist:
        if len(bashlist) == 1 and os.path.isdir(mask):
            for file in os.listdir(mask):
                path = os.path.join(mask, file)
                if os.path.isfile(path):
                    paths.append(path)
        else:
            for path in bashlist:
                if os.path.isfile(path):
                    paths.append(path)

    elif len(mask) == 4:
        paths.append(mask)

    return sorted(paths)


tsv_columns = ['R', 'Q', 'ID', 'RL', 'QL', 'TYPE', 'LALI',
               'RMSD', 'RTM', 'QTM','MATCH']
tsv_head    = '\t'.join(tsv_columns)
tsv_linefmt = '{R}\t{Q}\t{ID}\t{RL}\t{QL}\t{TYPE}\t{LALI}\t'\
              '{RMSD:6.2f}\t{RTM:6.5f}\t{QTM:6.5f}\t{MATCH}'

def rescode(mulind:'pd.MultiIndex'):
    return (mulind
            .to_series()
            .reset_index()
            .replace('?', '')
            .iloc[:,:5]
            .astype(str)
            .apply('.'.join, axis=1))


def tsv_line_data(artemis:'ARTEMIS', **data) -> 'tuple[str,str]':
    # SEQORD TOPIND
    if state == 0:      # FF
        if artemis.perm:
            ans = artemis.get_permutation()
            mch = ','.join(  rescode(artemis.ans2['rAli'])
                           + '='
                           + rescode(artemis.ans2['qAli']))

            data['TYPE']  = 'TOPIND'
            data['LALI']  = ans['p_aliLength']
            data['RMSD']  = ans['p_RMSD']
            data['RTM']   = ans['p_rTMscore']
            data['QTM']   = ans['p_qTMscore']
            data['MATCH'] = mch
        else:
            ans = artemis.get_alignment()
            ri, qi = hitFromAli(artemis.ans1['rAli'],
                                artemis.ans1['qAli'])
            mch = ','.join(  rescode(artemis.r.resi[ri])
                            + '='
                            + rescode(artemis.q.resi[qi]))

            data['TYPE']  = 'SEQORD'
            data['LALI']  = ans['aliLength']
            data['RMSD']  = ans['RMSD']
            data['RTM']   = ans['rTMscore']
            data['QTM']   = ans['qTMscore']
            data['MATCH'] = mch

        return data, 

    elif state == 1:    # FT
        ans = artemis.get_permutation()
        mch = ','.join(  rescode(artemis.ans2['rAli'])
                       + '='
                       + rescode(artemis.ans2['qAli']))

        data['TYPE']  = 'TOPIND'
        data['LALI']  = ans['p_aliLength']
        data['RMSD']  = ans['p_RMSD']
        data['RTM']   = ans['p_rTMscore']
        data['QTM']   = ans['p_qTMscore']
        data['MATCH'] = mch

        return data, 

    elif state == 2:    # TF
        ans = artemis.get_alignment()
        ri, qi = hitFromAli(artemis.ans1['rAli'],
                            artemis.ans1['qAli'])
        mch = ','.join(  rescode(artemis.r.resi[ri])
                       + '='
                       + rescode(artemis.q.resi[qi]))

        data['TYPE']  = 'SEQORD'
        data['LALI']  = ans['aliLength']
        data['RMSD']  = ans['RMSD']
        data['RTM']   = ans['rTMscore']
        data['QTM']   = ans['qTMscore']
        data['MATCH'] = mch

        return data,

    elif state == 3:    # TT
        data_2 = data.copy()

        ans = artemis.get_permutation()
        mch = ','.join(  rescode(artemis.ans2['rAli'])
                       + '='
                       + rescode(artemis.ans2['qAli']))

        data_2['TYPE']  = 'TOPIND'
        data_2['LALI']  = ans['p_aliLength']
        data_2['RMSD']  = ans['p_RMSD']
        data_2['RTM']   = ans['p_rTMscore']
        data_2['QTM']   = ans['p_qTMscore']
        data_2['MATCH'] = mch

        data_1 = data.copy()

        ans = artemis.get_alignment()
        ri, qi = hitFromAli(artemis.ans1['rAli'],
                            artemis.ans1['qAli'])
        mch = ','.join(  rescode(artemis.r.resi[ri])
                       + '='
                       + rescode(artemis.q.resi[qi]))

        data_1['TYPE']  = 'SEQORD'
        data_1['LALI']  = ans['aliLength']
        data_1['RMSD']  = ans['RMSD']
        data_1['RTM']   = ans['rTMscore']
        data_1['QTM']   = ans['qTMscore']
        data_1['MATCH'] = mch

        return data_1, data_2


def save(artemis:'ARTEMIS', saveto, ID):
    # SEQORD TOPIND
    if state == 0:      # FF
        if artemis.perm:
            artemis.save_topind(
                saveto      = saveto, 
                saveformat  = args.saveformat,
                saveres     = args.saveres,
                ID          = ID
            )

        else:
            artemis.save_seqord(
                saveto      = saveto, 
                saveformat  = args.saveformat,
                saveres     = args.saveres,
                ID          = ID
            )

    elif state == 1:    # FT
        artemis.save_topind(
            saveto      = saveto, 
            saveformat  = args.saveformat,
            saveres     = args.saveres,
            ID          = ID
        )

    elif state == 2:    # TF
        artemis.save_seqord(
            saveto      = saveto, 
            saveformat  = args.saveformat,
            saveres     = args.saveres,
            ID          = ID
        )

    elif state == 3:    # TT
        artemis.save_seqord(
            saveto      = saveto, 
            saveformat  = args.saveformat,
            saveres     = args.saveres,
            ID          = ID
        )
        artemis.save_topind(
            saveto      = saveto, 
            saveformat  = args.saveformat,
            saveres     = args.saveres,
            ID          = ID
        )


def rcut(artemis:'ARTEMIS') -> 'DataModel':

    r = artemis.r.copy()
    resi = r.resi

    # SEQORD TOPIND
    if state == 0:      # FF
        if artemis.perm: 
            rAli = artemis.ans2['rAli']
            resi = resi.difference(rAli,
                                   sort=False)
            seedi = r.seedi.difference(rAli,
                                       sort=False)

        else:
            ri, qi = hitFromAli(artemis.ans1['rAli'],
                                artemis.ans1['qAli'])
            rAli  = resi[ri]
            resi  = resi.difference(rAli,
                                    sort=False)
            seedi = r.seedi.difference(rAli,
                                       sort=False)

    elif state == 1:    # FT
        rAli = artemis.ans2['rAli']
        resi = resi.difference(rAli,
                               sort=False)
        seedi = r.seedi.difference(rAli,
                                   sort=False)

    elif state == 2:    # TF
        ri, qi = hitFromAli(artemis.ans1['rAli'],
                            artemis.ans1['qAli'])
        rAli  = resi[ri]
        resi  = resi.difference(rAli, sort=False)
        seedi = r.seedi.difference(rAli,
                                   sort=False)

    elif state == 3:    # TT
        ri, qi = hitFromAli(artemis.ans1['rAli'],
                            artemis.ans1['qAli'])
        rAli  = resi[ri]
        resi  = resi.difference(rAli,
                                sort=False)
        seedi = r.seedi.difference(rAli,
                                 sort=False)

        rAli = artemis.ans2['rAli']
        resi = resi.difference(rAli,
                               sort=False)
        seedi = seedi.difference(rAli,
                                 sort=False)



    r.resi = resi
    r.seedi = seedi

    s = r.resrepr1.loc[seedi]
    s = s[s.notna()].values
    r.s = s

    m = r.resrepr2.loc[resi]
    m = m[m.notna()]
    r.i = m.index

    if not m.empty:
        m = np.vstack(m.values)
        t = KDTree(m)
    else:
        m = np.array([])
        t = None
    r.m = m
    r.t = t

    return r


def getQTM(artemis:'ARTEMIS'):
    # SEQORD TOPIND
    if state == 0:      # FF
        if artemis.perm:
            return artemis.ans2['qTM']
        else:
            return artemis.ans1['qTM']
    elif state == 1:      # FT
        return artemis.ans2['qTM']
    elif state == 2:      # TF
        return artemis.ans1['qTM']
    elif state == 3:      # TT
        return max(artemis.ans1['qTM'], artemis.ans2['qTM']) 


if __name__ == '__main__':
    from src.argparse import argParse

    args = argParse(sys.argv[1:])

    if args.addhits > 0.0:
        args.tsv = True
        if args.addhits >= 1.0:
            args.addhits = int(args.addhits)


    if args.tsv:
        state = int('{:d}{:d}'.format(args.sequential, 
                                      args.permutation), base=2)
        print(tsv_head)

    if args.threads > 1:
        if 'fork' in mp.get_all_start_methods():
            mp.set_start_method('fork')

    rpaths = pathparser(args.r)
    qpaths = pathparser(args.q)

    if not rpaths:
        raise ValueError(
            'Argument r={} does not specify any local file or structure from the RCSB PDB.'
            .format(args.r)
        )

    if not qpaths:
        raise ValueError(
            'Argument q={} does not specify any local file or structure from the RCSB PDB.'
            .format(args.q)
        )


    for rpath in rpaths:
        try:
            rformat = ({'.cif': '.cif'}.get(os.path.splitext(rpath)[-1], '.pdb')
                    if args.rformat is None
                    else args.rformat)

            r = DataModel(
                path    = rpath,
                fmt     = rformat,
                res     = args.rres,
                resneg  = args.rresneg,
                seed    = args.rseed
            )

            if args.saveto:
                if len(rpaths) > 1:
                    _, file = os.path.split(rpath)
                    saveto = os.path.join(args.saveto, file)
                else:
                    saveto = args.saveto

        except Exception as e:
            if args.silent:
                print('Exception: r={}'.format(rpath), file=sys.stderr)
                tb = e.__traceback__
                print(''.join(traceback.format_tb(tb)), file=sys.stderr)
                continue
            else:
                raise e

        for qpath in qpaths:
            try:
                qformat = ({'.cif': '.cif'}.get(os.path.splitext(qpath)[-1], '.pdb')
                        if args.qformat is None
                        else args.qformat)

                q = DataModel(
                    path    = qpath,
                    fmt     = qformat,
                    res     = args.qres,
                    resneg  = args.qresneg,
                    seed    = args.qseed
                )

                artemis = ARTEMIS(
                    r, q,
                    matchrange  = args.matchrange,
                    threads     = args.threads,
                    toplargest  = args.toplargest,
                    shift       = args.shift,
                    stepdiv     = args.stepdiv,
                    superonly   = args.superonly,
                    notmopt     = args.notmopt,
                )
                artemis.run()

                if args.tsv:
                    data = {
                        'R' : rpath,
                        'RL': r.L,
                        'Q' : qpath,
                        'QL': q.L,
                        'ID': 0,
                    }

                    for data in tsv_line_data(artemis, **data):
                        print(tsv_linefmt.format(**data))

                    if args.saveto:
                        save(artemis, saveto, data['ID'])

                    if isinstance(args.addhits, int):
                        r_ = rcut(artemis)
                        data['RL'] = r_.L
                        for i in range(args.addhits):
                            if r_.L > q.L / 2:
                                data['ID'] += 1
                                artemis = ARTEMIS(
                                    r_, q,
                                    matchrange  = args.matchrange,
                                    threads     = args.threads,
                                    toplargest  = args.toplargest,
                                    shift       = args.shift,
                                    stepdiv     = args.stepdiv,
                                    superonly   = args.superonly,
                                    notmopt     = args.notmopt,
                                )
                                artemis.run()

                                for data in tsv_line_data(artemis, **data):
                                    print(tsv_linefmt.format(**data))
                                if args.saveto:
                                    save(artemis, saveto, data['ID'])

                                r_ = rcut(artemis)
                                data['RL'] = r_.L
                            else:
                                break

                    elif args.addhits > 0:
                        r_ = rcut(artemis)
                        data['RL'] = r_.L
                        QTM = getQTM(artemis)
                        while QTM > args.addhits and r_.L > q.L / 2:
                            data['ID'] += 1
                            artemis = ARTEMIS(
                                r_, q,
                                matchrange  = args.matchrange,
                                threads     = args.threads,
                                toplargest  = args.toplargest,
                                shift       = args.shift,
                                stepdiv     = args.stepdiv,
                                superonly   = args.superonly,
                                notmopt     = args.notmopt,
                            )
                            artemis.run()

                            for data in tsv_line_data(artemis, **data):
                                if data['QTM'] > args.addhits:
                                    print(tsv_linefmt.format(**data))
                                    if args.saveto:
                                        if data['TYPE'] == 'SEQORD':
                                            artemis.save_seqord(
                                                saveto      = saveto, 
                                                saveformat  = args.saveformat,
                                                saveres     = args.saveres,
                                                ID          = data['ID']
                                            )
                                        elif data['TYPE'] == 'TOPIND':
                                            artemis.save_topind(
                                                saveto      = saveto, 
                                                saveformat  = args.saveformat,
                                                saveres     = args.saveres,
                                                ID          = data['ID']
                                            )
                            r_ = rcut(artemis)
                            data['RL'] = r_.L
                            QTM = getQTM(artemis)

                else:
                    if args.saveto:
                        artemis.save(saveto      = saveto, 
                                    saveformat  = args.saveformat,
                                    saveres     = args.saveres,
                                    permutation = artemis.perm or
                                                  args.permutation)

                    print(artemis
                        .show(verbose     = args.verbose, 
                                permutation = args.permutation or
                                            artemis.perm,
                                dssr=args.dssr))

            except Exception as e:
                if args.silent:
                    print('Exception: r={}, q={}'.format(rpath, qpath), file=sys.stderr)
                    tb = e.__traceback__
                    print(''.join(traceback.format_tb(tb)), file=sys.stderr)
                    continue
                else:
                    raise e
