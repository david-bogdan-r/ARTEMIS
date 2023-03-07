import multiprocessing as mp
import os
import re
import sys
from itertools import product

import numpy as np
import pandas as pd
import tabulate
from scipy.spatial import KDTree
from src import NAR, PDB, artmsNW

BASEDIR = os.path.dirname(os.path.abspath(__file__))
BASEFMT = {'PDB', 'CIF', 'MMCIF'}
HELPARG = {'--H', '-H', '--h', '-h', '--help', '-help'}
COORD   = ['Cartn_x', 'Cartn_y', 'Cartn_z']
BASEREP = NAR.get_repr('5')



def help(path:'str'= BASEDIR + '/help.txt'):
    with open(path, 'r') as file:
        print(file.read())


def transform(x:'np.ndarray', y:'np.ndarray'):
    Ex = x.mean(axis=0)
    Ey = y.mean(axis=0)
    xx = x - Ex
    yy = y - Ey
    m  = np.dot(np.transpose(yy), xx)
    u, s, vh = np.linalg.svd(m)
    a = np.transpose(np.dot(np.transpose(vh), np.transpose(u)))
    if np.linalg.det(a) < 0:
        vh[2] = -vh[2]
        a = np.transpose(np.dot(np.transpose(vh), np.transpose(u)))
    b = Ex - np.dot(Ey, a)
    return a, b


def RMSD(x:'np.ndarray', y:'np.ndarray'):
    d = x - y
    return np.sqrt(np.sum(np.sum(np.multiply(d, d))) / len(d))


def MCS(we:'list'):
    v2w = {}
    v1w = {}
    for we in we:
        v1, v2, w = we
        if v1 not in v2w:
            v2w[v1] = v2, w
        else:
            if w < v2w[v1][1]:
                v2w[v1] = v2, w
        if v2 not in v1w:
            v1w[v2] = v1, w
        else:
            if w < v1w[v2][1]:
                v1w[v2] = v1, w
    vv = []
    for v2 in v1w.keys():
        v1, w = v1w[v2]
        if v2w[v1][0] == v2:
            vv.append((v1, v2))
    return vv


def resrep(res):
    base = res.atom_site['auth_comp_id'].iloc[0]
    if base not in BASEREP:
        return None
    else:
        atom_site = res.atom_site.copy()
        atom_site.set_index('auth_atom_id', inplace=True)
        
        mat = []
        for atom in BASEREP[base]:
            if atom == 'P' and atom not in atom_site.index:
                atom = "O5'"
            try:
                vec = atom_site.loc[atom.split(), COORD].values
                mat.append(vec.mean(axis=0))
            except:
                return None
        return np.array(mat)


def get_structure(key:'str') -> 'PDB.Structure': 
    global kwargs
    
    if key not in {'r', 'q'}:
        raise ValueError('{} not in {}'.format(key, {'r', 'q'}))
    
    path = str(kwargs.get(key))
    if path is None:
        raise ValueError('Key argument {}={} is {}'.format(key, path, None))
    
    fmt = kwargs.get(key + 'format')
    name, ext = os.path.splitext(os.path.basename(path))
    
    if fmt is None:
        fmt = ext[1:].upper()
        if fmt not in BASEFMT:
            fmt = 'PDB'
    else:
        fmt = fmt.upper()
    
    if fmt == 'MMCIF':
        fmt = 'CIF'

    if fmt == 'PDB':
        struct = PDB.read_pdb(path, name)
    else:
        struct = PDB.read_cif(path, name)
    
    struct.drop_altloc(keep='last')

    return struct


def parce(key:'str'):
    global kwargs
    
    struct = get_structure(key)
    
    res    = kwargs.get(key + 'res'   , '#1')
    resneg = kwargs.get(key + 'resneg', '')
    seed   = kwargs.get(key + 'seed'  , res)
    
    flg = seed != res
    
    neg = bool(resneg)
    substruct = struct.get_substruct([res, resneg][neg], neg)
    
    res   = []
    mat   = []
    seedi = []
    i     = -1    
    for r in substruct.resgen():
        x = resrep(r)
        if x is None:
            continue
        i += 1
        res.append(r)
        mat.append(x)
        if flg and r._res_mask(seed).any():
            seedi.append(i)
        else:
            seedi.append(i)
    
    mat   = np.array(mat)
    base  = [r.atom_site['auth_comp_id'].iloc[0] for r in res]
    seqex = ''.join([b if len(b) == 1 else '({})'.format(b) for b in base])
    seq   = (re.sub(r'\(\w+\)', 'm', seqex), )
    
    out = {
        'name'      : struct.name,
        'struct'    : struct,
        'substruct' : substruct,
        'res'       : res,
        'mat'       : mat,
        'seed'      : seedi,
        'base'      : base,
        'seq'       : seq,
        'seqex'     : seqex
    }
    return out


def gethit(seed):
    global rmat,  qmat
    global rtree, qmean
    
    i, j = seed
    x, y = rmat[i], qmat[j]
    a, b = transform(x, y)
    
    qtree = KDTree(np.dot(qmean, a) + b)
    d = rtree.sparse_distance_matrix(qtree,
                                     max_distance=matchrange,
                                     p=2,
                                     output_type='ndarray')
    mcs = MCS(d)
    if mcs:
        i, j = map(list, zip(*mcs))
        size = len(mcs)
    else:
        return None
    
    if (size < sizemin) or (size > sizemax):
        return None
        
    x, y = rmat[i], np.dot(qmat[j], a) + b
    rmsd = RMSD(np.vstack(x), np.vstack(y))
    if (rmsd < rmsdmin) or (rmsd > rmsdmax):
        return None
        
    rmsdsize = rmsd / size
    if (rmsdsize < rmsdsizemin) or (rmsdsize > rmsdsizemax):
        return None
    
    resrmsd = max(map(RMSD, x, y))
    if (resrmsd < resrmsdmin) or (resrmsd > resrmsdmax):
        return None
    
    return zip(i, j)


def fillmat(hits):
    global mat
    for hit in hits:
        if hit is None:
            continue
        for i, j in hit:
            mat[i, j] += 1


def asterisk(raseq, qaseq, mat, threshold, M) -> 'str':
    i, j = -1, -1
    a = ''
    for r, q in zip(raseq, qaseq):
        if r != '-':
            i += 1
        if q != '-':
            j += 1
        if mat[i, j] / M > threshold:
            a += '*'
        else:
            a += '.'
    return a



if __name__ == '__main__':
    args = sys.argv[1:]
    if any(arg in HELPARG for arg in args):
        help()
        exit()
    
    kwargs = dict(arg.split('=') for arg in args)
    
    threads = int(kwargs.get('threads', 1))
    if threads < 1 or threads > mp.cpu_count():
        threads = mp.cpu_count()
    if threads > 1:
        mp.set_start_method('fork')
    
    rmsdmin = float(kwargs.get('rmsdmin', 0.))
    rmsdmax = float(kwargs.get('rmsdmax', 1e10))
    
    sizemin = float(kwargs.get('sizemin', 1.))
    sizemax = float(kwargs.get('sizemax', 1e10))
    
    rmsdsizemin = float(kwargs.get('rmsdsizemin', 0.))
    rmsdsizemax = float(kwargs.get('rmsdsizemax', 1e10))
    
    resrmsdmin = float(kwargs.get('resrmsdmin', 0.))
    resrmsdmax = float(kwargs.get('resrmsdmax', 1e10))
    
    matchrange = float(kwargs.get('matchrange', 3.))
    
    threshold = float(kwargs.get('threshold', 0.5))
    if threshold < 0. or threshold > 1.:
        raise ValueError('threshold < 0.' if threshold < 0. 
                         else 'threshold > 1.')
    
    r = parce('r')
    q = parce('q')
    
    M     = min(len(r['seed']), len(q['seed']))
    seed  = product(r['seed'], q['seed'])
    
    rtree = KDTree(r['mat'].mean(axis=1))
    qmean = q['mat'].mean(axis=1)
    
    n, m  = len(r['res']), len(q['res'])
    mat   = np.zeros((n, m), dtype=int)
    rmat  = r['mat']
    qmat  = q['mat']
    if threads > 1:
        hits = mp.Pool(threads).map(gethit, seed)
        fillmat(hits)
    else:
        hits = list(map(gethit, seed))
        fillmat(hits)
    
    
    backtrack, score = artmsNW.GapGlobalAlign(r['seq'], q['seq'], sigma=0., eps=0., scoremat=mat)
    
    raseq = artmsNW.OutputGapAlign1(backtrack, 
                                    r['seq'], 
                                    len(r['seq'][0]), len(q['seq'][0]))[0]
    qaseq = artmsNW.OutputGapAlign2(backtrack, 
                                    q['seq'], 
                                    len(r['seq'][0]),
                                    len(q['seq'][0]))[0]
    
    i, j = -1, -1
    ast  = ''
    th   = threshold * M
    for ri, qi in zip(raseq, qaseq):
        if ri != '-':
            i += 1
        if qi != '-':
            j += 1
        if mat[i, j] > th:
            ast += '*'
        else:
            ast += '.'
    
    
    p = mat[mat != 0] / mat.sum()
    h = -(p * np.log2(p)).sum() / np.log2(p.size)
    
    qrespos  = [i for i, c in enumerate(qaseq) if c != '-']
    i = -1
    rresposrev = {}
    for ii, c in enumerate(raseq):
        if c != '-':
            i += 1
            rresposrev[ii] = i
        else:
            rresposrev[ii] = -1
    
    mats = mat.astype(str).T
    for i in range(len(qrespos)):
        j = rresposrev[qrespos[i]]
        mats[i,j] = '<' + mats[i,j] + '>'
    
    out = [
        ['IDENTIFIED REFERENCE SEQUENCE', r['seqex']],
        ['IDENTIFIED QUERY SEQUENCE', q['seqex']],
        ['REFERENCE', raseq],
        ['QUERY', qaseq],
        ['ASTERISK', ast],
        ['ASTERISKSCORE', ast.count('*') / min(m, n)],
        ['ENTROPY', h],
    ]
    
    if (r['name'], q['name']) in {('4kqy', '2gis'), ('2gis', '4kqy')}:
        if r['name'] == '4kqy':
            _raseq = 'GUUCUUAUCAAGAGA-AGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAAC'
            _qaseq = '-GGCUUAUCAAGA-GAGGUGGAGGGACUGGCCCGAUGAAACCC-GGCAACC----------AGAAAU----------GGUGCCAAUUCCUGCAGCG-GAAA--CGUUGAAAGAUGAGCCA'
        else:
            _qaseq = 'GUUCUUAUCAAGAGA-AGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAAC'
            _raseq = '-GGCUUAUCAAGA-GAGGUGGAGGGACUGGCCCGAUGAAACCC-GGCAACC----------AGAAAU----------GGUGCCAAUUCCUGCAGCG-GAAA--CGUUGAAAGAUGAGCCA'
        
        _qrespos = [i for i, c in enumerate(_qaseq) if c != '-']
        i = -1
        _rresposrev = {}
        for ii, c in enumerate(_raseq):
            if c != '-':
                i += 1
                _rresposrev[ii] = i
            else:
                _rresposrev[ii] = -1
        
        lattice = ''
        for i in range(len(qrespos)):
            if i == 0:
                lattice += '-' * (_qrespos[i])
            else:
                lattice += '-' * (_qrespos[i] - _qrespos[i-1] - 1)
            if rresposrev[qrespos[i]] == _rresposrev[_qrespos[i]]:
                lattice += '#'
            else:
                lattice += '.'
        out += [
            ['RREF', _raseq],
            ['QREF', _qaseq],
            ['LATTICE', lattice],
            ['LATTICESCORE', lattice.count('#') / (lattice.count('#') + lattice.count('.'))]
        ]
        
        for i in range(len(_qrespos)):
            j = _rresposrev[_qrespos[i]]
            mats[i,j] = '[' + mats[i,j] + ']'
        out = tabulate.tabulate(out, tablefmt='plain')
    mats = pd.DataFrame(mats, columns=r['base'])
    mats.index = q['base']
    
    mats.to_csv(sys.stdout, sep='\t')
    print(out)