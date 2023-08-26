import numpy as np

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
