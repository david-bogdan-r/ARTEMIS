from typing import Callable

import numpy as np
from scipy.spatial.distance import cdist

from src.Hit import closest, hitFromAli
from src.Kabsch import transform
from src.NW import globalAlign

def Impose(
    rm:'np.ndarray', qm:'np.ndarray',
    rL:'int', qL:'int',
    rd2:'float', qd2:'float'
):

    '''
    Return (rot, tran), (rTM, qTM), RMSD of optimal superposition
    between rm and qm matrices
    '''

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

    return opt(0, N) # (rot, tran), (rTM, qTM), RMSD


def Align(
    ri:'np.ndarray', qi:'np.ndarray',
    shift:'float',
    rc:'np.ndarray', qc:'np.ndarray',
    rm:'np.ndarray', qm:'np.ndarray',
    rseq:'str', qseq:'str',
    impose:'Callable'
):

    mask:'np.ndarray' = np.array([rc[i] is not None for i in ri])\
        | np.array([qc[i] is not None for i in qi])
    mask = np.where(mask)[0]

    rx  = np.vstack(rc[ri[mask]])
    qx  = np.vstack(qc[qi[mask]])

    a, b = transform(rx, qx)
    dist = cdist(rm, np.dot(qm, a) + b)

    ri, qi = np.where(dist < 8)
    h = closest(dict(zip(zip(ri, qi), dist[ri, qi])))
    ri, qi = h[:, h[0].argsort()]
    s1 = impose(rm[ri], qm[qi])
    ans1 = (*s1, (ri, qi))

    scoremat = -dist
    scoremat -= scoremat.min()
    scoremat -= min(scoremat[ri, qi]) - shift

    rAli, qAli = globalAlign(rseq, qseq, scoremat)
    ri, qi = hitFromAli(rAli, qAli)
    s2 = impose(rm[ri], qm[qi])
    ans2 = (*s2, (rAli, qAli))

    return ans1, ans2
