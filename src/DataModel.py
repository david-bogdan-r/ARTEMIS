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
    res:'pd.Series'

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
            [b if len(b) == 1 else 'M' 
             for b in self.i.get_level_values(2)]
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

    def copy(self):

        new = DataModel(
            name=self.name, fmt=self.fmt, 
            atom_site=self.atom_site.copy(),
            i=self.i.copy(), 
            c=self.c.copy(), 
            m=self.m.copy(), 
            t=self.t,
            res=self.res,
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
    resNeg:'bool',
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

    r     = resRepr.apply(atom_site)
    rmcbi = r.reset_index()[MCBI].values

    res   = r.loc[getResSpec(rmcbi, resSpec , resNeg)]
    notna = res.notna()

    if not notna.any():
        raise Exception(
            '{}. No residue was found in (res={}, resneg={})'.format(name, resSpec, resNeg)
        )

    atom_site_c = atom_site.set_index(MCBI).loc[res[notna].index]
    atom_site_c.set_index('auth_atom_id', inplace=True)

    i = res[notna].index
    c = res[notna].values
    m:'np.ndarray' = atom_site_c.loc["C3'", CRDN].values # type: ignore
    if len(m) < 3:
        raise Exception(
            '{}. Not enough residue for alignment (res={}, resneg={})'.format(name, resSpec, resNeg)
        )
    t = KDTree(m)

    if seedSpec == resSpec and resNeg:
        loc  = getResSpec(rmcbi, seedSpec, resNeg)
    else:
        loc = getResSpec(rmcbi, seedSpec, False)

    if not loc:
        raise Exception(
            '{}. Seed residues is empty (seed={})'.format(name, seedSpec)
        )
    seed = r.loc[loc]


    saveres = getResSpec(rmcbi, saveresSpec, False) 

    dm = DataModel(
        name=name, fmt=fmt, atom_site=atom_site,
        i=i, c=c, m=m, t=t, res=res, seed=seed, saveres=saveres # type: ignore
    )

    return dm
