import json

import pandas as pd
import numpy as np

MCBI     = [
    'pdbx_PDB_model_num',
    'auth_asym_id',
    'auth_comp_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code'
]
CRDN     = ['Cartn_x', 'Cartn_y', 'Cartn_z']
DROPCOL  = [
    'pdbx_PDB_model_num',
    'auth_asym_id',
    'auth_comp_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code',
    'auth_atom_id'
]
KEEP    = 'last'


def resrepr(atom_site:'pd.DataFrame', **kwargs:'list[str]') -> 'pd.Series':

    def repr(res_atom_site:'pd.DataFrame'):
        base  = res_atom_site['auth_comp_id'].iloc[0]

        if base not in resr:
            return None

        currep = resr[base]
        mat = []

        for atom in currep:
            inter = atom.intersection(res_atom_site.index)

            if inter.empty:
                return None

            coord = res_atom_site.loc[inter, CRDN].values
            if len(inter) > 1:
                coord = coord.mean(axis=0)

            mat.append(coord)

        return np.vstack(mat)

    atom_site = (atom_site
                 .drop_duplicates(DROPCOL, keep=KEEP)
                 .set_index('auth_atom_id'))

    resr:'dict[str, list[pd.Index]]' = {}

    for k, v in kwargs.items():
        i = [pd.Index(a.split()) for a in v]
        resr[k] = i

    rep = (atom_site
           .groupby(MCBI, sort=False)
           .apply(repr)) # type: ignore

    return rep

def load_resrepr(*paths):
    reprs   = []
    passed  = {}

    for path in paths:
        if path not in passed:
            with open(path, 'r') as file:
                passed[path] = json.load(file)
        reprs.append(passed[path])

    return reprs
