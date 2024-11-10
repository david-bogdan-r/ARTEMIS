import os
import warnings
from copy import deepcopy
from io import IOBase
from string import ascii_letters, digits

import numpy as np
import pandas as pd


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
        '{occupancy:>6}',
        '{B_iso_or_equiv:>6}          ',
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
    'auth_seq_id',
    'pdbx_PDB_ins_code'
]
CRDN = ['Cartn_x', 'Cartn_y', 'Cartn_z']

URL = 'https://files.rcsb.org/view/{PDB_ID}{fmt}'

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

    atom_site = pd.DataFrame(a).astype(
        {
            'id':int,
            'Cartn_x': float,
            'Cartn_y': float,
            'Cartn_z': float,
            'auth_seq_id': int,
            'pdbx_PDB_model_num': int
        }
    )

    if 'auth_atom_id' in atom_site.columns:
        atom_site['auth_atom_id'] = atom_site['auth_atom_id'].str.replace('*', "'")
        atom_site['auth_atom_id'] = atom_site['auth_atom_id'].replace(
            to_replace={
                'O1P': 'OP1',
                'O2P': 'OP2'
            }
        )

    return atom_site

def as_cif(text:'str') -> 'pd.DataFrame':

    s = text.find('loop_\n_atom_site.') + len('loop_\n')
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
    atom_site = pd.DataFrame(a, columns=columns).astype(
        {
            'id':int,
            'Cartn_x': float,
            'Cartn_y': float,
            'Cartn_z': float,
            'auth_seq_id': int,
            'pdbx_PDB_model_num': int
        }
    )

    for col in atom_site.columns:
        if col.startswith('label'):
            new_col = col.replace('label', 'auth')
            if new_col not in atom_site.columns:
                atom_site[new_col] = atom_site[col]

    for col in ['label_atom_id', 'auth_atom_id']:
        if col in atom_site.columns:
            m = atom_site[col].str.startswith('"', "'")
            atom_site.loc[m, col] = atom_site[col][m].str[1:-1]

    if 'auth_atom_id' in atom_site.columns:
        atom_site['auth_atom_id'] = atom_site['auth_atom_id'].str.replace('*', "'")
        atom_site['auth_atom_id'] = atom_site['auth_atom_id'].replace(
            to_replace={
                'O1P': 'OP1',
                'O2P': 'OP2'
            }
        )
    if 'label_atom_id' in atom_site.columns:
        atom_site['label_atom_id'] = atom_site['label_atom_id'].str.replace('*', "'")
        atom_site['label_atom_id'] = atom_site['label_atom_id'].replace(
            to_replace={
                'O1P': 'OP1',
                'O2P': 'OP2'
            }
        )

    return atom_site


def read_pdb(path:'str|IOBase') -> 'pd.DataFrame':
    atom_site = as_pdb(read(path))

    if atom_site['Cartn_x'].dtype != float:
        if isinstance(path, str):
            raise Exception('It looks like the file={} read incorrectly.\n'
                            'Check that the specified file format is correct.'
                            .format(path))


    return atom_site


def read_cif(path:'str|IOBase') -> 'pd.DataFrame':
    atom_site = as_cif(read(path))

    if atom_site.empty:
        raise Exception('It looks like the file={} read incorrectly.\n'
                        'Check that the specified file format is correct.'
                        .format(path))

    return atom_site


class BaseModel:

    def __init__(self, path:'str', fmt:'str') -> 'None':

        self.path = path
        self.fmt  = fmt

        name, _ = os.path.splitext(path)
        name = name.split(os.sep)[-1]

        if fmt == '.cif':
            if os.path.isfile(path):
                atom_site = read_cif(path)
            else:
                if len(name) == 4:
                    import requests
                    url = URL.format(PDB_ID=name.upper(), fmt='.cif')
                    response = requests.get(url)
                    if response.status_code == 200:
                        atom_site = as_cif(response.text)
                    else:
                        raise Exception(
                            'Could not obtain {name} structure in {fmt} format from RCSB PDB'.format(
                                name=name, fmt=fmt.lstrip('.').upper()
                            )
                        )
                else:
                    raise ValueError(f"path={path} it's not a local file or PDB ID [4 letter code]")
        else:
            if os.path.isfile(path):
                atom_site = read_pdb(path)
            else:
                if len(name) == 4:
                    import requests
                    url = URL.format(PDB_ID=name.upper(), fmt='.pdb')
                    response = requests.get(url)
                    if response.status_code == 200:
                        atom_site = as_pdb(response.text)
                    else:
                        raise Exception(
                            'Could not obtain {name} structure in {fmt} format from RCSB PDB'.format(
                                name=name, fmt=fmt.lstrip('.').upper()
                            )
                        )
                else:
                    raise ValueError(f"path={path} it's not a local file or PDB ID [4 letter code]")

        self.path = path
        self.name = name
        self.fmt  = fmt
        self.atom_site = atom_site


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
        return deepcopy(self)

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
            atom_site.replace('', '?', inplace=True)

            prea = 'auth_'
            prel = 'label_'
            for col in ['seq_id', 'comp_id', 'asym_id', 'atom_id']:
                cola = prea + col
                coll = prel + col 
                atom_site[coll] = atom_site[cola]

            label = atom_site['label_asym_id'].unique()
            label = dict(zip(label, range(1, len(label) + 1)))

            with pd.option_context("future.no_silent_downcasting", True):
                atom_site['label_entity_id'] = (atom_site['label_asym_id']
                                                .replace(to_replace=label))

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

            label_seq_id = atom_site[[
                'pdbx_PDB_model_num',
                'auth_asym_id', 
                'auth_comp_id',
                'auth_seq_id',
                'pdbx_PDB_ins_code'
            ]].astype(str).apply(lambda x: '.'.join(x), axis=1)
            res_id = label_seq_id.unique()
            replace = dict(zip(res_id, range(1, len(res_id)+1)))
            with pd.option_context("future.no_silent_downcasting", True):
                label_seq_id.replace(replace, inplace=True)
            atom_site['label_seq_id'] = label_seq_id

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

            cols = [
                'occupancy',
                'B_iso_or_equiv',
                'type_symbol',
                'pdbx_formal_charge',
                'pdbx_PDB_model_num'
            ]
            for col in cols:
                if col not in atom_site.columns:
                    atom_site[col] = ''

        return atom_site


    def get_cif_text(self) -> 'str':

        atom_site = self.get_cif_atom_site()
        atom_site = add_atom_id_quote(atom_site)
        atom_site[CRDN] = atom_site[CRDN].round(3)

        text = 'data_{}\n'.format(str(self).upper())
        text += '# \n'
        text += 'loop_\n'

        for col in atom_site.columns:
            text += '_atom_site.{} \n'.format(col)

        colen = atom_site.apply(lambda x: max(x.astype(str).map(len))) + 1
        line_format = ''
        for col, l in colen.items():
            line_format += '{' + f'{col}:<{l}' + '}'
        text += '\n'.join(map(lambda x: line_format.format(**x[1]),
                              atom_site.iterrows()))
        text += '\n# \n'

        return text

    def get_pdb_text(self) -> 'str':

        atom_site = self.get_pdb_atom_site()
        if atom_site.empty:
            return ''

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


        if atom_site['occupancy'].dtype == float:
            atom_site['occupancy'] = atom_site['occupancy'].apply(
                lambda x: '{:.2f}'.format(x)
            )
        if atom_site['B_iso_or_equiv'].dtype == float:
            atom_site['B_iso_or_equiv'] = atom_site['B_iso_or_equiv'].apply(
                lambda x: '{:.2f}'.format(x)
            )
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
                text += PDBFORMAT['TER'].format(**a0) # type: ignore
                text += 'ENDMDL\n'

                m0 = a['pdbx_PDB_model_num']
                c0 = a['auth_asym_id']
                ii = 1
                text += PDBFORMAT['MODEL'].format(m0)

            elif a['auth_asym_id'] != c0:
                c0  = a['auth_asym_id']
                ii += 1
                a0['id'] = ii % 100_000
                text += PDBFORMAT['TER'].format(**a0) # type: ignore

            ii += 1
            a['id'] = ii % 100_000
            text += PDBFORMAT['ATOM'].format(**a) # type: ignore
            a0 = a

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
    rename = {c: '' for c in code if len(c) > 2}

    if rename:

        lbl1 = ascii_letters + digits
        lbl = (x + y for x in ' ' + lbl1 for y in lbl1)
        i = 0

        for c in rename.keys():
            try:
                l = next(lbl).lstrip()
                while l in code:
                    l = next(lbl).lstrip()
                rename[c] = l
                i += 1

            except IndexError:
                raise NameError('Unable to rename chains '\
                                'to a two-letter code')

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

    elif isinstance(type_symbol, str) and type_symbol and atom_id.startswith(type_symbol):
        suff = atom_id[len(type_symbol):]
        return '{:>2}'.format(type_symbol) + '{:<2}'.format(suff)

    else:
        if len(atom_id) == 1:
            return atom_id+' '*2
        if len(atom_id) == 2:
            return atom_id+' '
        return atom_id

def specParse(spec:'str') -> 'list':

    def parse(spec:'str') -> 'tuple':

        model:'int|None|str' = ''
        chain:'str|None'     = ''
        base :'str|None'     = ''
        start:'int|None|str' = ''
        end  :'int|None|str' = ''
        ins_code = ''

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
                model = 1

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

        if not base:
            base = None

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

                start = int(dg)
                ins_code = ch

            else:
                start = int(start)
        else:
            start = None


        return model, chain, base, start, ins_code, end

    return [parse(s) for s in spec.split()]

def getResSpec(atom_site:'pd.DataFrame', spec:'str', neg:'bool'=False) -> 'list':

    res = atom_site[MCBI].drop_duplicates().values
    n   = len(res)
    msk = np.array([False] * n)

    res = res.T

    for m, c, b, s, ins, e in specParse(spec):

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
                cur_msk &= (res[4] == ins)

        msk |= cur_msk
    res = res.T

    if neg:
        msk ^= True

    return [tuple(r) for r in res[msk]]
