import os
from pathlib import Path

import numpy as np
import pandas as pd

from string import ascii_letters, digits

BASEDIR = Path(__file__).parent.as_posix()
COORD   = ['Cartn_x', 'Cartn_y', 'Cartn_z']
PDBFMT  = {
    'ATOM'  :   '{group_PDB:<6}{id:>5} {auth_atom_id:<4}'\
                '{label_alt_id:1}{auth_comp_id:>3}{auth_asym_id:>2}'\
                '{auth_seq_id:>4}{pdbx_PDB_ins_code:1}   '\
                '{Cartn_x:>8.3f}{Cartn_y:>8.3f}{Cartn_z:>8.3f}'\
                '{occupancy:>6.2f}{B_iso_or_equiv:>6.2f}          '\
                '{type_symbol:>2}{pdbx_formal_charge:>2}\n',
    'TER'   :   'TER   {id:>5}      {auth_comp_id:>3}{auth_asym_id:>2}'\
                '{auth_seq_id:>4}\n',
    'MODEL' :   'MODEL     {:>4}\n',
    'REMARK':   'REMARK 250 CHAIN RENAMING {} -> {}\n',
    'ENDMDL':   'ENDMDL\n'
}

pd.to_numeric.__defaults__ = 'ignore', None

class Structure:

    def __init__(self, name:'str', atom_site:'pd.DataFrame'):
        self.name = name
        self.atom_site = atom_site

    def __str__(self):
        return self.name

    def __repr__(self):
        return '<{} Structure>'.format(self)

    def to_cif(self, path:'str'):

        folder, file = os.path.split(path)
        if folder:
            os.makedirs(folder, exist_ok=True)

        txt = ''.join(('data_{}\n# \nloop_\n'.format(self),
                       '\n'.join(['_atom_site.{}'.format(col)
                                  for col in self.atom_site.columns]) + '\n',
                       self.atom_site.to_string(header=False, index=False),
                       '\n# \n'))

        with open(path, 'w') as file:
            file.write(txt)

    def to_pdb(self, path:'str'):

        folder, file = os.path.split(path)
        if folder:
            os.makedirs(folder, exist_ok=True)

        atom_site = self.atom_site.copy()
        atom_site.replace('.', '', inplace=True)
        atom_site.replace('?', '', inplace=True)

        with open(path, 'w') as file:
            chains = atom_site['auth_asym_id'].astype(str).unique()
            rename = {c: '' for c in chains if len(c) > 1}
            if rename:
                i = 0
                n = ascii_letters + digits
                for c in rename.keys():
                    try:
                        while n[i] in rename:
                            i += 1
                        rename[c] = n[i]
                        i += 1
                    except IndexError:
                        raise NameError('Unable to rename chains'\
                                        'to a one-letter code')
                for kv in rename.items():
                    file.write(PDBFMT['REMARK'].format(*kv))
                atom_site.replace({'auth_asym_id':rename}, inplace=True)

            m = atom_site.iloc[0]['pdbx_PDB_model_num']
            c = atom_site.iloc[0]['auth_asym_id']
            dc = 0
            dm = 0
            a = {}

            file.write(PDBFMT['MODEL'].format(m))
            for atom in atom_site.to_dict('index').values():
                if atom['pdbx_PDB_model_num'] != m:
                    file.write(PDBFMT['ENDMDL'])
                    m = atom['pdbx_PDB_model_num']
                    file.write(PDBFMT['MODEL'].format(m))

                elif atom['auth_asym_id'] != c:
                    c = atom['auth_asym_id']
                    a['id'] += 1
                    file.write(PDBFMT['TER'].format(**a))
                    dc += 1

                atom['id'] = (atom['id'] + dc + dm) % 100_000
                a = atom
                file.write(PDBFMT['ATOM'].format(**atom))

class ResRepr:

    def __init__(self, key:'str') -> 'None':
        self.key  = key
        self.repr = ResRepr.__read_repr(key)

    def __repr__(self):
        return '<ResRepr({})>'.format(self.key)

    def __getitem__(self, key:'str'):
        return self.repr[key]

    @staticmethod
    def __read_repr(key:'str') -> 'dict[str, str]':

        path = BASEDIR + '/ResRepr/{}.json'.format(key)
        with open(path, 'r') as file:
            d = eval(file.read())

        return d

    def __call__(self, atom_site:'pd.DataFrame') -> 'np.ndarray':

        if atom_site.index.name != 'auth_atom_id':
            atom_site = atom_site.set_index('auth_atom_id')

        resrepr = self.repr
        base = atom_site['auth_comp_id'].iloc[0]
        if base not in resrepr.keys():
            raise KeyError('Base {} not in {}/resrepr/{}.json.'
                            .format(base, self.key))

        m = []
        for a in resrepr[base]:
            if ' ' in a:
                al = a.split()
                an = len(al)
            else:
                if a == 'P' and a not in atom_site.index:
                    a = "O5'"
                al = [a]
                an = 1

            try:
                v = atom_site.loc[al, COORD]
                if an > 1:
                    v = v.mean(axis=0)
                m.append(v)
            except KeyError:
                aa = [ai for ai in a if ai not in atom_site.index]
                raise KeyError('atom_id {} not in atom_site.auth_atom_id.'
                               .format(aa))

        m = np.vstack(m)

        return m

def read_cif(path:'str', name:'str'='') -> 'Structure':

    if not name:
        name = os.path.splitext(os.path.basename(path))[0]

    with open(path, 'r') as file:
        t = file.read()

    col = t.find('_atom_site.')
    e = t.find('#', col) - 1
    t = t[col:e].split('\n')

    col = []
    i = 0
    l = t[i]
    while l.startswith('_'):
        col.append(l.split('.')[1].rstrip())
        i += 1
        l = t[i]

    t = map(str.split, t[i:])
    atom_site = pd.DataFrame(t, columns=col)
    atom_site = atom_site.apply(pd.to_numeric)

    col = 'label_atom_id'
    q = atom_site[col].str.startswith('"', '"')
    if q.any():
        atom_site.loc[q, col] = atom_site[q][col].str[1:-1]

    col = 'auth_atom_id'
    if col in atom_site:
        q = atom_site[col].str.startswith('"', '"')
        if q.any():
            atom_site.loc[q, col] = atom_site[q][col].str[1:-1]

    pa = 'auth_'
    pl = 'label_'
    for col in ['seq_id', 'comp_id', 'asym_id', 'atom_id']:
        ca = pa + col
        cl = pl + col
        if ca not in atom_site and cl in atom_site:
            atom_site[ca] = atom_site[cl]

    return Structure(name, atom_site)

def read_pdb(path:'str', name:'str'='') -> 'Structure':

    if not name:
        name = os.path.splitext(os.path.basename(path))[0]

    p = {'ATOM  ', 'HETATM'}
    m = '1'
    a = []
    with open(path, 'r') as file:
        for l in file:
            lp = l[:6]
            if lp in p:
                la = [
                    l[0 : 6].strip(),
                    l[6 :11].strip(),
                    l[12:16].strip(),
                    l[16:17].strip(),
                    l[17:20].strip(),
                    l[20:22].strip(),
                    l[22:26].strip(),
                    l[26:27].strip(),
                    l[30:38].strip(),
                    l[38:46].strip(),
                    l[46:54].strip(),
                    l[54:60].strip(),
                    l[60:66].strip(),
                    l[76:78].strip(),
                    l[78:80].strip()
                ]
                atom_id = la[2]
                if atom_id.startswith('"') or atom_id.startswith("'"):
                    la[2] = atom_id[1:-1]
                la.append(m)
                a.append(la)
            elif lp == 'MODEL ':
                m = l.split()[1]

    atom_site = pd.DataFrame(
        a,
        columns=[
            'group_PDB',
            'id',
            'auth_atom_id',
            'label_alt_id',
            'auth_comp_id',
            'auth_asym_id',
            'auth_seq_id',
            'pdbx_PDB_ins_code',
            'Cartn_x',
            'Cartn_y',
            'Cartn_z',
            'occupancy',
            'B_iso_or_equiv',
            'type_symbol',
            'pdbx_formal_charge',
            'pdbx_PDB_model_num'
        ]
    )
    atom_site = atom_site.apply(pd.to_numeric)

    atom_site['pdbx_PDB_ins_code'].fillna('?', inplace=True)
    atom_site['pdbx_formal_charge'].fillna('?', inplace=True)
    atom_site['label_alt_id'].fillna('.', inplace=True)

    pa = 'auth_'
    pl = 'label_'
    for col in ['seq_id', 'comp_id', 'asym_id', 'atom_id']:
        ca = pa + col
        cl = pl + col
        atom_site[cl] = atom_site[ca]

    entity_id = atom_site['label_asym_id'].unique()
    entity_id = dict(zip(entity_id, range(1, len(entity_id) + 1)))
    atom_site['label_entity_id'] = (atom_site['label_asym_id']
                                    .replace(entity_id))

    atom_site = atom_site[
            [
                'group_PDB',
                'id',
                'type_symbol',
                'label_atom_id',
                'label_alt_id',
                'label_comp_id',
                'label_asym_id',
                'label_entity_id',
                'label_seq_id',
                'pdbx_PDB_ins_code',
                'Cartn_x',
                'Cartn_y',
                'Cartn_z',
                'occupancy',
                'B_iso_or_equiv',
                'pdbx_formal_charge',
                'auth_seq_id',
                'auth_comp_id',
                'auth_asym_id',
                'auth_atom_id',
                'pdbx_PDB_model_num'
            ]
        ]

    return Structure(name, atom_site)
