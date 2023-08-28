from dataclasses import dataclass
from io import IOBase

import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from src.PDBio import CRDN, FORMAT, MCBI, Model, getResSpec, read_cif, read_pdb
from src.ResRepr import ResRepr


@dataclass
class DataModel(Model):

    i:'pd.Index'        # PDBio.MCBI ~ DSSR code
    c:'np.ndarray'      # Residue representation
    m:'np.ndarray'      # C3' coordinates
    t:'KDTree'          # KDTree of C3' coordinates

    seed:'pd.Series'
    saveres:'list'

    def __str__(self) -> 'str':
        return self.name

    def __repr__(self) -> 'str':
        return '<{} DataModel>'.format(self)

    def __len__(self) -> 'int':
        return len(self.i)

    def get_seq(self) -> 'str':

        seq = ''.join(
            [b if len(b) == 1 else 'm' 
             for b in self.i.get_level_values(2)]
        ).lower()

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

    def copy(self):

        new = DataModel(
            name=self.name, fmt=self.fmt, 
            atom_site=self.atom_site.copy(),
            i=self.i.copy(), 
            c=self.c.copy(), 
            m=self.m.copy(), 
            t=self.t,
            seed=self.seed,
            saveres=self.saveres
        )

        return new


    L   = property(__len__)
    seq = property(get_seq)
    d0  = property(get_d0)


def dataModel(
    name:'str', fmt:'str', path:'str|IOBase',
    resSpec:'str', seedSpec:'str', saveresSpec:'str', 
    resneg:'bool',
    resRepr:'ResRepr',
) -> 'DataModel':

    fmt = fmt.strip('.').upper()
    fmt = FORMAT.get(fmt, '.pdb')

    atom_site = (read_pdb, read_cif)[fmt == '.cif'](path)

    atom_site.drop_duplicates(
        [
            'pdbx_PDB_model_num',
            'auth_asym_id',
            'auth_comp_id',
            'auth_seq_id',
            'auth_atom_id'
        ],
        keep='last',
        inplace=True
    )

    atom_site_ = atom_site[atom_site['auth_comp_id'].isin(resRepr.repr.keys())]
    r      = resRepr.apply(atom_site_)
    rarray = r.reset_index()[MCBI].values

    res     = r.loc[getResSpec(rarray, resSpec , resneg)]

    atom_site_ = atom_site.set_index(MCBI).loc[res.index]
    atom_site_.set_index('auth_atom_id', inplace=True)

    i = res.index
    c = res.values
    m = atom_site_.loc["C3'", CRDN].values # type: ignore
    t = KDTree(m)

    if seedSpec == resSpec and resneg:
        seed = r.loc[
            getResSpec(rarray, seedSpec, resneg)
        ]
    else:
        seed = r.loc[
            getResSpec(rarray, seedSpec, False)
        ]

    saveres = getResSpec(rarray, saveresSpec, False) 

    dm = DataModel(
        name=name, fmt=fmt, atom_site=atom_site,
        i=i, c=c, m=m, t=t, seed=seed, saveres=saveres # type: ignore
    )

    return dm

