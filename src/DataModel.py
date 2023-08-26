from dataclasses import dataclass
from io import IOBase

import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from src.PDBio import CRDN, FORMAT, MCBI, Model, read_cif, read_pdb
from src.ResRepr import ResRepr


@dataclass
class DataModel(Model):

    i:'pd.Index'        # PDBio.MCBI ~ DSSR code
    c:'np.ndarray'      # Residue representation
    m:'np.ndarray'      # C3' coordinates
    t:'KDTree'          # KDTree of C3' coordinates

    def __str__(self) -> 'str':
        return self.name

    def __repr__(self) -> 'str':
        return '<{} DataModel>'.format(self)

    def __len__(self) -> 'int':
        return len(self.i)

    def get_seq(self) -> 'str':
        return ''.join(self.i.get_level_values(2)).lower()

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


def dataModel(
    name:'str', fmt:'str', path:'str|IOBase', resRepr:'ResRepr',
) -> 'DataModel':

    fmt = fmt.strip('.').upper()
    fmt = FORMAT.get(fmt, '.pdb')

    atom_site = (read_pdb, read_cif)[fmt == '.cif'](path)

    r = resRepr.apply(atom_site)
    r.dropna(inplace=True)

    atom_site_ = atom_site.set_index(MCBI).loc[r.index]
    atom_site_.set_index('auth_atom_id', inplace=True)

    i = r.index
    c = r.values
    m = atom_site_.loc["C3'", CRDN].values # type: ignore
    t = KDTree(m)

    dm = DataModel(
        name=name, fmt=fmt, atom_site=atom_site, 
        i=i, c=c, m=m, t=t  # type: ignore
    )

    return dm
