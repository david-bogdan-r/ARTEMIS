import numpy as np
from scipy.spatial import KDTree

from src.Kabsch import transform


def getHit(
    rx:'np.ndarray', qx:'np.ndarray',
    rt:'KDTree',     qm:'np.ndarray',
    matchrange
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


def closest(hit:'dict') -> 'np.ndarray':

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


def closestHit(
    rx:'np.ndarray', qx:'np.ndarray',
    rt:'KDTree', qm:'np.ndarray',
    matchrange
) -> 'np.ndarray':

    return closest(getHit(rx, qx, rt, qm, matchrange))

def mutuallyClosestHit(
    rx:'np.ndarray', qx:'np.ndarray',
    rt:'KDTree', qm:'np.ndarray',
    matchrange
) -> 'np.ndarray':
    
    return mutuallyClosest(getHit(rx, qx, rt, qm, matchrange))


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
