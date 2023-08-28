import os
from dataclasses import dataclass
from io import IOBase
from string import ascii_letters, digits

import numpy as np
import pandas as pd

pd.to_numeric.__defaults__ = 'ignore', None


FORMAT = {
    'PDB'  : '.pdb',
    'CIF'  : '.cif',
    'MMCIF': '.cif'
}

PDBFORMAT = {
    'ATOM'  : ''.join([
        '{group_PDB:<6}',
        '{id:>5} ',
        '{auth_atom_id:>4}',
        '{label_alt_id:1}',
        '{auth_comp_id:>3}',
        '{auth_asym_id:>2}',
        '{auth_seq_id:>4}',
        '{pdbx_PDB_ins_code:1}   ',
        '{Cartn_x:>8.3f}',
        '{Cartn_y:>8.3f}',
        '{Cartn_z:>8.3f}',
        '{occupancy:>6.2f}',
        '{B_iso_or_equiv:>6.2f}          ',
        '{type_symbol:>2}',
        '{pdbx_formal_charge:>2}\n'
    ]),
    'TER'   : ''.join([
        'TER   {id:>5}      ',
        '{auth_comp_id:>3}',
        '{auth_asym_id:>2}',
        '{auth_seq_id:>4}\n',
    ]),
    'MODEL' : 'MODEL     {:>4}\n',
    'REMARK': 'REMARK 250 CHAIN RENAMING {} -> {}\n',
}

MCBI = [
    'pdbx_PDB_model_num',
    'auth_asym_id',
    'auth_comp_id',
    'auth_seq_id'
]
CRDN = ['Cartn_x', 'Cartn_y', 'Cartn_z']


def read(path:'str|IOBase') -> 'str':

    if isinstance(path, str) and os.path.isfile(path):
        with open(path, 'r') as file:
            text = file.read()

    elif isinstance(path, IOBase):
        text = path.read()

    else:
        raise ValueError('path != path to file | readable object')

    return text


def as_pdb(text:'str') -> 'pd.DataFrame':

    columns = [
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

    k = {'ATOM  ', 'HETATM'}
    m = '1'
    a = []

    for l in text.split('\n'):

        lk = l[:6]
        if lk in k:
            al = [
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

            atom_id = al[2]
            if atom_id.startswith('"') or atom_id.startswith("'"):
                al[2] = atom_id[1:-1]

            al.append(m)
            a.append(dict(zip(columns, al)))

        elif lk == 'MODEL ':
            m = l.split()[1]

    atom_site = pd.DataFrame(a).apply(pd.to_numeric)
    atom_site['auth_asym_id'] = atom_site['auth_asym_id'].astype(str)

    return atom_site

def as_cif(text:'str') -> 'pd.DataFrame':

    s = text.find('_atom_site.')
    e = text.find('#', s) - 1
    t = text[s:e].split('\n')

    columns = []
    i = 0
    l = t[i]
    while l.startswith('_'):
        columns.append(l.split('.')[1].rstrip())
        i += 1
        l = t[i]

    a = map(str.split, t[i:])
    atom_site = pd.DataFrame(a, columns=columns)
    atom_site = atom_site.apply(pd.to_numeric)
    atom_site['auth_asym_id'] = atom_site['auth_asym_id'].astype(str)

    for col in ['label_atom_id', 'auth_atom_id']:
        if col in atom_site.columns:
            m = atom_site[col].str.startswith('"', "'")
            atom_site.loc[m, col] = atom_site[col][m].str[1:-1]

    return atom_site


def read_pdb(path:'str|IOBase') -> 'pd.DataFrame':
    return as_pdb(read(path))

def read_cif(path:'str|IOBase') -> 'pd.DataFrame':
    return as_cif(read(path))


@dataclass
class Model:

    name:'str'
    fmt :'str'
    atom_site:'pd.DataFrame'


    def __str__(self) -> 'str':
        return self.name

    def __repr__(self) -> 'str':
        return '<{}{} Model>'.format(self.name, self.fmt)


    def to_cif(self, path:'str')->'None':

        mkdir(path)

        text = self.get_cif_text()
        with open(path, 'w') as file:
            file.write(text)

    def to_pdb(self, path:'str')->'None':

        mkdir(path)

        text = self.get_pdb_text()
        with open(path, 'w') as file:
            file.write(text)


    def copy(self):

        new = Model(
            self.name,
            self.fmt,
            self.atom_site.copy()
        )

        return new


    def get_cif_atom_site(self) -> 'pd.DataFrame':

        atom_site = self.atom_site.copy()

        if self.fmt == '.pdb':
            atom_site.fillna(
                {
                    'label_alt_id'      : '.',
                    'pdbx_PDB_ins_code' : '?',
                    'pdbx_formal_charge': '?',
                    'type_symbol'       : atom_site['auth_atom_id'].str[0],
                    'occupancy'         : 1.,
                    'B_iso_or_equiv'    : 0.,
                },
                inplace=True
            )

            prea = 'auth_'
            prel = 'label_'
            for col in ['seq_id', 'comp_id', 'asym_id', 'atom_id']:
                cola = prea + col
                coll = prel + col 
                atom_site[coll] = atom_site[cola]

            label = atom_site['label_asym_id'].unique()
            label = dict(zip(label, range(1, len(label) + 1)))

            atom_site['label_entity_id'] = (atom_site['label_asym_id']
                                            .replace(label))

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

        return atom_site

    def get_pdb_atom_site(self) -> 'pd.DataFrame':

        atom_site = self.atom_site.copy()

        if self.fmt == '.cif':

            prea = 'auth_'
            prel = 'label_'

            for col in ['seq_id', 'comp_id', 'asym_id', 'atom_id']:
                cola = prea + col

                if cola not in atom_site.columns:
                    coll = prel + col
                    atom_site[cola] = atom_site[coll]

        return atom_site


    def get_cif_text(self) -> 'str':

        atom_site = self.get_cif_atom_site()
        atom_site = add_atom_id_quote(atom_site)

        text = 'data_{}\n'.format(str(self).upper())
        text += '# \n'
        text += 'loop_\n'

        for col in atom_site.columns:
            text += '_atom_site.{} \n'.format(col)

        for al in atom_site.astype(str).values:
            text += ' '.join(al) + '\n'

        text += '# \n'

        return text

    def get_pdb_text(self) -> 'str':

        atom_site = self.get_pdb_atom_site()
        chain_names = chain_rename(atom_site)
        atom_site['auth_atom_id'] = atom_site.apply(atom_id_PDBformat, axis=1)

        atom_site.replace(
            {
                'label_alt_id'      : {'.':''},
                'pdbx_PDB_ins_code' : {'?':''},
                'pdbx_formal_charge': {'?':''},
                'auth_asym_id': chain_names,
            },
            inplace=True
        )

        atom_site.fillna(
            {
                'type_symbol'       : atom_site['auth_atom_id'].str[0],
                'occupancy'         : 1.,
                'B_iso_or_equiv'    : 0.,
                'pdbx_formal_charge': '',
                'label_alt_id'      : '',
                'pdbx_PDB_ins_code' : ''
            },
            inplace=True
        )

        text = ''
        if chain_names:
            for i in chain_names.items():
                text += PDBFORMAT['REMARK'].format(*i)

        records = atom_site.to_dict('records')
        a0 = records[0]
        m0 = a0['pdbx_PDB_model_num']
        text += PDBFORMAT['MODEL'].format(m0)
        c0 = a0['auth_asym_id']
        ii = 0
        a  = {} # for type checking
        for a in records:
            if a['pdbx_PDB_model_num'] != m0:

                ii += 1
                a0['id'] = ii % 100_000
                text += Model.PDBFORMAT['TER'].format(**a0) # type: ignore
                text += 'ENDMDL\n'

                m0 = a['pdbx_PDB_model_num']
                c0 = a['auth_asym_id']
                ii = 1
                text += PDBFORMAT['MODEL'].format(m0)

            elif a['auth_asym_id'] != c0:
                c0  = a['auth_asym_id']
                ii += 1
                a0['id'] = ii % 100_000
                text += Model.PDBFORMAT['TER'].format(**a0) # type: ignore

            ii += 1
            a['id'] = ii % 100_000
            text += PDBFORMAT['ATOM'].format(**a) # type: ignore

        ii += 1
        a['id'] = ii % 100_000
        text += PDBFORMAT['TER'].format(**a) # type: ignore
        text += 'ENDMDL\n'
        text += 'END\n'

        return text


    def get_coord(self)->'np.ndarray':
        return self.atom_site[CRDN].values

    def set_coord(self, coord:'np.ndarray')->'None':
        self.atom_site.loc[:, CRDN] = coord

    coord = property(get_coord, set_coord)


def mkdir(path:'str') -> 'None':

    folder = os.path.split(path)[0]

    if folder:
        os.makedirs(folder, exist_ok=True)

def chain_rename(atom_site:'pd.DataFrame') -> 'dict':

    code = atom_site['auth_asym_id'].astype(str).unique()
    rename = {c: '' for c in code if len(c) > 1}

    if rename:

        lbl = ascii_letters + digits
        i = 0

        for c in rename.keys():
            try:
                while lbl[i] in code:
                    i += 1
                rename[c] = lbl[i]
                i += 1

            except IndexError:
                raise NameError('Unable to rename chains '\
                                'to a one-letter code')

    return rename

def add_atom_id_quote(atom_site:'pd.DataFrame')->'pd.DataFrame':

    for col in ['label_atom_id', 'auth_atom_id']:
        if col in atom_site.columns:

            m = atom_site[col].str.startswith('"', "'")

            if (True ^ m).all():
                atom_site[col] = '"' + atom_site[col] + '"'

            else:
                atom_site.loc[m, col] = atom_site[col][m].str[1:-1]
                atom_site[col] = '"' + atom_site[col] + '"'

    return atom_site

def atom_id_PDBformat(atom:'pd.Series') -> 'str':

    type_symbol = atom['type_symbol']
    atom_id     = atom['auth_atom_id']

    if len(atom_id) == 4:
        return atom_id

    elif isinstance(type_symbol, str) and type_symbol:
        suff = atom_id[len(type_symbol):]
        return '{:>2}'.format(type_symbol) + '{:<2}'.format(suff)

    else:
        return atom_id

def specParse(spec:'str') -> 'list':

    def parse(spec:'str') -> 'tuple':

        model:'int|None|str' = ''
        chain:'str|None'     = ''
        base :'str|None'     = ''
        start:'int|None|str' = ''
        end  :'int|None|str' = ''

        resid:'str'          = ''

        model_is_None = True

        state = -1
        for c in spec:
            if c == '#':
                state = 0
                model_is_None = False
                continue
            elif c == '/':
                state = 1
                continue
            elif c == ':':
                state = 2
                continue

            if state == 0:
                model += c
            elif state == 1:
                chain += c
            elif state == 2:
                resid += c

        if model:
            model = int(model)
        else:
            if model_is_None:
                model = None

        if not chain:
            chain = None

        state = 0
        for c in resid:
            if c == '_':
                state += 1
                continue

            if state == 0:
                base += c

            elif state == 1:
                start += c

            elif state == 2:
                end += c

        if not end:
            end = None
        else:
            end = int(end)

        if start:
            if not start.isdigit():
                ch = ''
                dg = ''

                for c in start:
                    if c.isdigit():
                        dg += c
                    else:
                        ch += c

                start = dg

                if not base and end is None:
                    base  = ch

            else:
                start = int(start)
        else:
            start = None

        if not base:
            base = None


        return model, chain, base, start, end

    return [parse(s) for s in spec.split()]

def getResSpec(res:'np.ndarray', spec:'str', neg:'bool'=False) -> 'list':

    n   = len(res)
    msk = np.array([False] * n)

    res = res.T
    for m, c, b, s, e in specParse(spec):

        cur_msk = np.array([True] * n)

        if m is not None:
            if isinstance(m, int):
                cur_msk &= res[0] == m

        if c is not None:
            cur_msk &= res[1] == c

        if b is not None:
            cur_msk &= res[2] == b

        if s is not None:
            if e is not None:
                cur_msk &= (s <= res[3]) & (res[3] <= e)
            else:
                cur_msk &= res[3] == s

        msk |= cur_msk
    res = res.T

    if neg:
        msk ^= True

    return [tuple(r) for r in res[msk]]