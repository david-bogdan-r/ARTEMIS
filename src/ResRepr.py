import json
import os

import pandas as pd
import numpy as np

from src.PDBio import Model, MCBI, CRDN

PATH = os.path.dirname(__file__) + '/ResRepr/{}.json'

class ResRepr:

    def __init__(self, resrepr:'str') -> 'None':

        self.name = resrepr

        self.repr:'dict[str, list[pd.Index]]' = {}
        self.atom:'dict[str, pd.Index]'       = {}

        with open(PATH.format(resrepr), 'r') as file:
            d:'dict[str, list[str]]' = json.load(file)

        for k, v in d.items():

            i = [pd.Index(a.split()) for a in v]

            self.repr[k] = i

            a = pd.Index([])
            for ii in i:
                a = a.append(ii)

            self.atom[k] = a

    def __str__(self) -> 'str':
        return self.name

    def __repr__(self) -> 'str':
        return '<{} ResRepr>'.format(self)

    def resRepr(self, atom_site:'pd.DataFrame') -> 'np.ndarray|None':

        base    = atom_site['auth_comp_id'].iloc[0]
        resrepr = self.repr.get(base)

        if not resrepr:
            return None

        mat = []

        for atom in resrepr:
            inter = atom.intersection(atom_site.index)

            if inter.empty:
                return None

            coord = atom_site.loc[inter, CRDN].values

            if len(inter) > 1:
                coord = coord.mean(axis=0)

            mat.append(coord)

        return np.vstack(mat)

    def apply(self, model:'Model|pd.DataFrame') -> 'pd.Series':

        if isinstance(model, Model):
            atom_site = model.atom_site.set_index('auth_atom_id')

        else:
            atom_site = model.set_index('auth_atom_id')

        res = (atom_site
              .groupby(MCBI, sort=False)
              .apply(self.resRepr)) # type: ignore

        return res
