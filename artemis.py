import itertools
import json
import multiprocessing as mp
import os
import sys
from copy import deepcopy
from datetime import datetime as dt
from functools import partial
from heapq import nlargest
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist


from src.argparse import argParse
from src.Kabsch import transform
from src.NW import globalAlign
from src.PDBio import BaseModel, getResSpec
from src.resrepr import resrepr

BASEDIR  = os.path.dirname(__file__)
SEEDPOOL = 100_000
RESREPR1 = BASEDIR + '/src/resrepr/artemis.json'
RESREPR2 = "C3'"
MCBI     = [
    'pdbx_PDB_model_num',
    'auth_asym_id',
    'auth_comp_id',
    'auth_seq_id'
]
CRDN     = ['Cartn_x', 'Cartn_y', 'Cartn_z']

INDEX = '{head}{config}{alignment}{distance}{permutation}{time}'

HEAD ='''
 ********************************************************************
 * ARTEMIS (Version 20230828)                                       *
 * using ARTEM to Infer Sequence alignment                          *
 * Reference: TODO                                                  *
 * Please email comments and suggestions to dbohdan@iimcb.gov.pl    *
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

saveto={saveto}
saveformat={saveformat}
saveres={saveres}

matchrange={matchrange}
nlargest={nlargest}
shift={shift}
stepdiv={stepdiv}

threads={threads}
_____________________________________________________________________
'''

ALIGNMENT = '''
Name of structure r: {rName}:{rChain}
Name of structure q: {qName}:{qChain} (to be superimposed onto structure r)
Length of structure r: {rLength} residues
Length of structure q: {qLength} residues

Aligned length= {aliLength}, RMSD= {RMSD:6.2f}, Seq_ID=n_identical/n_aligned= {Seq_ID:4.3f}
TM-score= {rTMscore:6.5f} (normalized by length of structure r: L={rLength}, d0={r_d0:.2f})
TM-score= {qTMscore:6.5f} (normalized by length of structure q: L={qLength}, d0={q_d0:.2f})

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
{rAlignment}
{distances}
{qAlignment}
'''

DISTANCE = [
    '\n\n'
    'Distance table:\n'
    '{rName:<{L1}}  dist   {qName:<{L2}}\n',
    '{rRes:<{L1}}  {d:<5.2f}  {qRes:<{L2}}\n'
]

PERMUTATION = '''
_____________________________________________________________________
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

class DataModel(BaseModel):

    def __init__(self, path:'str', fmt:'str',
                repr1:'dict', repr2:'dict',
                res:'str', resneg:'str|None', seed:'str|None') -> 'None':

        super().__init__(path, fmt)

        self.repr1 = repr1
        self.repr2 = repr2

        self.res  = res
        self.resneg = resneg

        atom_site = self.atom_site

        resrepr1 = resrepr(atom_site, **repr1)
        resrepr2 = resrepr(atom_site, **repr2)

        self.resrepr1 = resrepr1
        self.resrepr2 = resrepr2

        atom_site = atom_site[atom_site['auth_comp_id']
                                .isin( list(repr1.keys()) \
                                      +list(repr2.keys()))]

        if resneg is None:
            resi = getResSpec(atom_site, res)
            erkey = 'res'
            erval   = res
            if seed is None:
                self.seed = self.res
            else:
                self.seed = seed
        else:
            resi = getResSpec(atom_site, resneg, True)
            erkey = 'resneg'
            erval = resneg
            if seed is None:
                self.seed = self.resneg
            else:
                self.seed = seed

        resi = pd.MultiIndex.from_tuples(resi, names=MCBI)
        self.resi  = resi

        if resi.empty:
            raise ValueError(
                'Empty residue set path={}:{}={}.'
                .format(path, erkey, erval)
            )

        if seed is None:
            seedi = resi
        else:
            seedi = getResSpec(atom_site, seed)
        seedi = pd.MultiIndex.from_tuples(seedi, names=MCBI)
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
                 r:'str', rformat:'str',
                 rres:'str', rresneg:'str|None', rseed:'str|None',
                 q:'str', qformat:'str',
                 qres:'str', qresneg:'str|None', qseed:'str|None',
                 matchrange:'float|None', threads:'int|None',
                 nlargest:'int|None', shift:'float|None', stepdiv:'int|None',
                 saveto:'str|None', saveformat:'str', saveres:'str', 
                 verbose:'bool', permutation:'bool',
                 repr1:'dict', repr2:'dict',
                 ) -> 'None':

        self.tic = dt.now()

        self.saveto     = saveto
        self.saveformat = saveformat
        self.saveres    = saveres

        self.verbose = verbose
        self.permutation = permutation

        if matchrange is None:
            matchrange = 3.5

        self.matchrange = matchrange

        if threads is None:
            threads = mp.cpu_count()

        self.threads = threads
        self.pool    = mp.Pool(threads)

        self.r = DataModel(
            path=r, fmt=rformat,
            repr1=repr1, repr2=repr2,
            res=rres, resneg=rresneg, seed=rseed
        )

        self.q = DataModel(
            path=q, fmt=qformat,
            repr1=repr1, repr2=repr2,
            res=qres, resneg=qresneg, seed=qseed
        )

        if self.q.L >= 500:
            if nlargest is None:
                nlargest = 2 * threads

            if shift is None:
                shift = 20

            if stepdiv is None:
                stepdiv = 100

        else:
            if nlargest is None:
                nlargest = self.q.L

            if shift is None:
                shift = 3

            if stepdiv is None:
                stepdiv = 0

        self.nlargest:'int' = nlargest # type: ignore
        self.shift    = shift
        self.stepdiv  = stepdiv

        self.hit = partial(
            self.Hit,
            rt=self.r.t, qm=self.q.m, matchrange=self.matchrange # type: ignore
        )

        self.impose = partial(
            self.Impose, 
            rL=self.r.L, qL=self.q.L,
            rd2=self.r.d0 ** 2, qd2=self.q.d0 ** 2
        )

        self.align = partial(
            self.Align,
            r=self.r, q=self.q,
            shift=self.shift,
            impose=self.impose
        )

    @staticmethod
    def Hit(
        rx:'np.ndarray', qx:'np.ndarray',
        rt:'KDTree', qm:'np.ndarray',
        matchrange:'float'
    ) -> 'dict':

        a, b = transform(rx, qx)
        qt = KDTree(np.dot(qm, a) + b)

        h = rt.sparse_distance_matrix(
            qt,
            max_distance=matchrange,
            p=2,
            output_type='dict'
        )

        return h

    @staticmethod
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

    @staticmethod
    def Closest(hit:'dict') -> 'np.ndarray':
        r = set()
        q = set()
        k = []

        for i, j in sorted(hit, key=hit.get): # type: ignore

            if i in r or j in q:
                continue

            else:
                r.add(i)
                q.add(j)
                k.append([i, j])

        h = np.array(k).T

        return h

    @staticmethod
    def mutuallyClosest(hit:'dict') -> 'np.ndarray':
        r = {}
        q = {}

        for i, j in hit:

            if i in r:
                r[i] += 1
            else:
                r[i] = 1

            if j in q:
                q[j] += 1
            else:
                q[j] = 1

        h = np.array([
            [i, j] for i, j in hit
            if r[i] == 1 and q[j] == 1
        ]).T

        return h

    @staticmethod
    def Impose(
        rm:'np.ndarray', qm:'np.ndarray',
        rL:'int', qL:'int',
        rd2:'float', qd2:'float'
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
        o = opt(0, N)

        ans = {
            'a'  : o[0][0],
            'b'  : o[0][1],
            'rTM': o[1][0],
            'qTM': o[1][1],
            'RMSD': o[2]
        }

        return ans

    def insertNaNBase(self, rAli, qAli, distances):

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
                distances.insert(j, ' ')
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
                distances.insert(j, ' ')
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
                    c = 'M'
                rAli += c

        qAli = ''
        for i in qali:
            if i == -1:
                qAli += '-'
            else:
                c = q.resi[i][2]
                if len(c) > 1:
                    c = 'M'
                qAli += c

        distances = ''.join(distances)

        return rAli, qAli, distances

    @staticmethod
    def Align(
        ri:'np.ndarray', qi:'np.ndarray',
        r:'DataModel', q:'DataModel',
        shift:'float',
        impose:'Callable'
    ):

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

        ri, qi = np.where(dist < 8)

        h = ARTEMIS.Closest(dict(zip(zip(ri, qi), dist[ri, qi])))
        ri, qi = h[:, h[0].argsort()]

        ans1 = impose(rm[ri], qm[qi])
        ans1['rAli'] = r.i[ri]
        ans1['qAli'] = q.i[qi]

        scoremat = -dist
        scoremat -= scoremat.min()
        scoremat -= min(scoremat[ri, qi]) - shift

        rAli, qAli = globalAlign(r.seq, q.seq, scoremat)
        h = ARTEMIS.hitFromAli(rAli, qAli)

        if h.size < 6:
            a = np.eye(3)
            b = np.zeros(3)

            ans2 = {
                'a': a,
                'b': b,
                'rTM': -1,
                'qTM': -1,
                'RMSD': -1,
                'rAli': rAli,
                'qAli': qAli,
            }

        else:
            ri, qi = h

            ans2 = impose(rm[ri], qm[qi])
            ans2['rAli'] = rAli
            ans2['qAli'] = qAli

        return ans1, ans2

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

        seed, count = self.get_seed()

        i = 0
        h = []

        while i + SEEDPOOL < count:
            step = [next(seed) for _ in range(SEEDPOOL)]
            h = nlargest(
                self.nlargest,
                h + self.get_hit(step),
                key=len
            )
            i += SEEDPOOL

        h = nlargest(
            self.nlargest,
            h + self.get_hit(seed),
            key=len
        )

        h = map(self.mutuallyClosest, h)
        h = [hh for hh in h if hh.size >= 6]

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

        pAns = max(
            ali,
            key=lambda x: x[0]['rTM'] + x[0]['qTM']
        )[0]

        Ans = max(
            ali,
            key=lambda x: x[1]['rTM'] + x[1]['qTM']
        )[1]

        if (pAns['qTM'] - Ans['qTM']) / abs(Ans['qTM']) >= 0.1:
            self.permutation = True

        self.pAns = pAns
        self.Ans  = Ans

    def get_config(self):

        r = self.r
        q = self.q

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
            'nlargest'  : self.nlargest,
            'shift'     : self.shift,
            'stepdiv'   : self.stepdiv,

            'saveto'    : self.saveto if self.saveto  else '',
            'saveformat': self.saveformat,
            'saveres'   : self.saveres,

            'threads'   : self.threads,

            'verbose'   : self.verbose,
            'permutation':self.permutation
        }

        return config

    def get_alignment(self) -> 'dict':

        r = self.r
        q = self.q

        Ans  = self.Ans

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

        rAli = Ans['rAli']
        qAli = Ans['qAli']

        aliLength = sum(i != '-' and j != '-' for i, j in zip(rAli, qAli))

        if aliLength:
            n_identical = sum(c1 == c2 for c1, c2 in zip(rAli, qAli))
            n_aligned   = aliLength
            Seq_ID = n_identical / n_aligned

            a = Ans['a']
            b = Ans['b']

            h = self.hitFromAli(rAli, qAli)
            rm = r.m[h[0]]
            qm = np.dot(q.m[h[1]], a) + b
            d = np.sqrt(((rm - qm) ** 2).sum(axis=1))

            distances = []
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
                        distances.append(':')
                    else:
                        distances.append('.')
                else:
                    distances.append(' ')

        else:
            Seq_ID = 0
            distances = [' '] * len(rAli)

        rAli, qAli, distances = (self
                                 .insertNaNBase(Ans['rAli'], 
                                                Ans['qAli'], 
                                                distances))

        alignment = {
            'rName'     : r.path,
            'qName'     : q.path,
            'rChain'    : rChain,
            'qChain'    : qChain,
            'rLength'   : r.L,
            'qLength'   : q.L,
            'aliLength' : aliLength,
            'RMSD'      : Ans['RMSD'],
            'Seq_ID'    : Seq_ID,
            'rTMscore'  : Ans['rTM'],
            'r_d0'      : r.d0,
            'qTMscore'  : Ans['qTM'],
            'q_d0'      : q.d0,
            'rAlignment': rAli,
            'qAlignment': qAli,
            'distances' : distances
        }

        return alignment

    def get_permutation(self) -> 'dict':

        Ans = self.pAns

        p_aliLength = len(Ans['rAli'])
        p_rTMscore  = Ans['rTM']
        p_qTMscore  = Ans['qTM']
        p_RMSD      = Ans['RMSD']

        p_n_identical = 0
        for rc, qc in zip(Ans['rAli'], Ans['qAli']):
            if rc[2] == qc[2]:
                p_n_identical += 1
        p_Seq_ID = p_n_identical / p_aliLength

        permutation = {
            'p_aliLength'   : p_aliLength,
            'p_rTMscore'    : p_rTMscore,
            'p_qTMscore'    : p_qTMscore,
            'p_RMSD'        : p_RMSD,
            'p_Seq_ID'      : p_Seq_ID
        }

        return permutation

    def get_distance(self):

        r = self.r
        q = self.q

        title = {
            'rName': r.name,
            'qName': q.name
        }

        Ans  = self.Ans
        rAli = Ans['rAli']
        qAli = Ans['qAli']

        aliLength = sum(i != '-' and j != '-' for i, j in zip(rAli, qAli))

        L1 = len(r.name)
        L2 = len(q.name)

        rows = []
        if aliLength:
            a = Ans['a']
            b = Ans['b']

            h = self.hitFromAli(rAli, qAli)
            rm = r.m[h[0]]
            qm = np.dot(q.m[h[1]], a) + b
            d = np.sqrt(((rm - qm) ** 2).sum(axis=1))

            for rb, qb, dd in zip(r.i[h[0]], q.i[h[1]], d):
                rRes = '.'.join(map(str, rb))
                qRes = '.'.join(map(str, qb))

                if len(rRes) > L1:
                    L1 = len(rRes)
                if len(qRes) > L2:
                    L2 = len(qRes)

                entity = {
                    'rRes'  : rRes,
                    'd'     : dd,
                    'qRes'  : qRes,
                }
                rows.append(entity)

            for entity in rows:
                entity['L1'] = L1
                entity['L2'] = L2

        title['L1'] = L1 # type: ignore
        title['L2'] = L2 # type: ignore
            
        return title, rows

    def get_time(self) -> 'dict':

        toc = dt.now()
        sec = (toc - self.tic).total_seconds()

        time = {
            'total_time': sec
        }

        return time

    def get_ans(self) -> 'str':

        index = {
            'head': HEAD,
            'alignment': ALIGNMENT.format(**self.get_alignment()),
            'time': TIME.format(**self.get_time()),
        }

        if self.permutation:
            index['permutation'] = PERMUTATION.format(**self.get_permutation())
        else:
            index['permutation'] = ''

        if self.verbose:
            index['config']   = CONFIG.format(**self.get_config())
            t, r = self.get_distance()
            s = DISTANCE[0].format(**t)
            rowpat = DISTANCE[1]
            for rr in r:
                s += rowpat.format(**rr)
            index['distance'] = s
        else:
            index['config']   = ''
            index['distance'] = ''

        return INDEX.format(**index)

    def save(self):

        if self.saveto is None:
            return

        saveto  = self.saveto.strip(os.sep)
        saveres = self.saveres
        saveformat = self.saveformat

        r = self.r
        q = self.q

        saveres = getResSpec(q.atom_site, saveres)

        Ans  = self.Ans
        pAns = self.pAns

        if os.path.isdir(saveto): # type: ignore
            files = set(os.listdir(saveto))
        else:
            files = set()

        rAli = Ans['rAli']
        qAli = Ans['qAli']

        h  = self.hitFromAli(rAli, qAli)

        if h.size:
            a = Ans['a']
            b = Ans['b']

            qq = deepcopy(q)
            qq.coord = np.dot(qq.coord, a) + b
            qq.atom_site = (qq.atom_site
                            .set_index(MCBI).loc[saveres]
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

            text = '{}\tdist\t{}\n'.format(r.name, q.name)
            for i, (ii, jj) in enumerate(zip(ri, qi)):
                text += "{}\t{:.2f}\t{}\n".format('.'.join(list(map(str, ii))),
                                                 d[i],
                                                 '.'.join(list(map(str, jj))))

            with open('{}/{}'.format(saveto, fname), 'w') as file:
                file.write(text)

        if self.permutation:
            a = pAns['a']
            b = pAns['b']

            qq = deepcopy(q)
            qq.coord = np.dot(qq.coord, a) + b
            qq.atom_site = (qq.atom_site
                            .set_index(MCBI).loc[saveres]
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

            ri = pAns['rAli'].to_list()
            qi = pAns['qAli'].to_list()

            rm = np.vstack(r.resrepr2.loc[ri])
            qm = np.vstack(q.resrepr2.loc[qi])
            qm = np.dot(qm, a) + b

            d = np.sqrt(((rm - qm) ** 2).sum(axis=1))

            fname = '{}_to_{}_p.tsv'.format(q.name, r.name)
            i = 0
            while fname in files:
                fname = '{}_to_{}_p_({}).tsv'.format(q.name, r.name, i)
                i += 1

            text = '{}\tdist\t{}\n'.format(self.r, self.q)
            for i, (ii, jj) in enumerate(zip(ri, qi)):
                text += "{}\t{:.2f}\t{}\n".format('.'.join(list(map(str, ii))),
                                                 d[i],
                                                 '.'.join(list(map(str, jj))))

            with open('{}/{}'.format(saveto, fname), 'w') as file:
                file.write(text)



if __name__ == '__main__':

    args = argParse(sys.argv[1:])

    if args.threads > 1:
        if 'fork' in mp.get_all_start_methods():
            mp.set_start_method('fork')

    with open(RESREPR1, 'r') as file:
        repr1 = json.load(file)
        repr2 = {k:RESREPR2.split() for k in repr1}

    artemis = ARTEMIS(**dict(args._get_kwargs()), 
                      repr1=repr1, repr2=repr2)

    artemis.run()
    artemis.save()
    print(artemis.get_ans())
