import numpy as np
from scipy.spatial.transform import Rotation

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

def wtransform(x:'np.ndarray', y:'np.ndarray', w:'np.ndarray'):

    Ex = x.mean(axis=0)
    Ey = y.mean(axis=0)

    xx = x - Ex
    yy = y - Ey

    a = Rotation.align_vectors(xx, yy, w)[0].as_matrix().T
    b = Ex - np.dot(Ey, a)

    return a, b
