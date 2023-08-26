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

    def __resRepr(self, atom_site:'pd.DataFrame') -> 'np.ndarray|None':

        base    = atom_site.name[2]
        resrepr = self.repr.get(base)

        if not resrepr:
            return None

        diff = self.atom[base].difference(atom_site.index)

        if not diff.empty:
            if ((len(diff) == 1)
                and (diff[0] == 'P')
                and ("O5'" in atom_site.index)):

                atom_site.loc['P'] = atom_site.loc[b"O5'"]

            else:
                return None

        mat = []

        for atom in resrepr:
            coord = atom_site.loc[atom, CRDN].values

            if len(atom) > 1:
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
              .apply(self.__resRepr)) # type: ignore

        return res
