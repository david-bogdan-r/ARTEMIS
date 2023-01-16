import os 
import functools

import pandas as pd
import numpy  as np


pd.set_option('mode.chained_assignment', None)
pd.to_numeric.__defaults__ = 'ignore', None


class structure(pd.DataFrame):
    counter = 0
    
    def __init__(self, *args, **kwargs) -> 'None':
        structure.counter += 1
        
        if not 'name' in kwargs:
            name = 'structure.{}'.format(structure.counter)
        else:
            name = kwargs['name']
            del kwargs['name']
        
        if not 'fmt' in kwargs:
            fmt = None
        else:
            fmt = kwargs['fmt']
            del kwargs['fmt']
        
        super().__init__(*args, **kwargs)
        
        self.name = name
        self.fmt  = fmt
    
    def __str__(self) -> 'str':
        return self.name
    
    def __repr__(self) -> str:
        return '<{} structure>'.format(self)
    
    def rename(self, name:'str') -> 'None':
        self.name = name
    
    def drop_altloc(self, keep:'str' = 'last'):
        super().drop_duplicates(
            [
                'pdbx_PDB_model_num', 
                'auth_asym_id', 
                'auth_comp_id',
                'auth_seq_id',
                'pdbx_PDB_ins_code',
                'auth_atom_id'
            ],
            keep=keep,
            inplace=True
        )
    
    
    def deprecate(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            other = func(self, *args, **kwargs)
            other = structure(other)
            other.name = self.name
            other.fmt  = self.fmt
            return other
        return wrapper
    
    @deprecate
    def apply(self, *args, **kwargs):
        return super().apply(*args, **kwargs)
    
    @deprecate
    def fillna(self, *args, **kwargs):
        return super().fillna(*args, **kwargs)
    
    @deprecate
    def copy(self, *args, **kwargs):
        return super().copy(*args, **kwargs)
    
    
    @property
    def coord(self):
        return self[['Cartn_x', 'Cartn_y', 'Cartn_z']]
    
    @coord.setter
    def coord(self, new_coord):
        self[['Cartn_x', 'Cartn_y', 'Cartn_z']] = new_coord
    
    @property
    def dssr_columns(self):
        return [
            'pdbx_PDB_model_num',
            'auth_asym_id',
            'auth_comp_id',
            'auth_seq_id',
            'pdbx_PDB_ins_code',
        ]
    
    @property
    def dssr(self):
        return self[self.dssr_columns].replace('?', '')
    
    
    def as_cif(self) -> 'structure':
        if self.fmt == 'CIF':
            return self.copy()
        
        tab = pd.DataFrame(self)
        tab['pdbx_PDB_ins_code'].replace('', '?', inplace=True)
        tab['pdbx_formal_charge'].replace('', '?', inplace=True)
        tab['label_alt_id'].replace('', '.', inplace=True)
        
        tab['label_atom_id']   = tab['auth_atom_id']
        tab['label_comp_id']   = tab['auth_comp_id']
        tab['label_asym_id']   = tab['auth_asym_id']
        tab['label_seq_id']    = tab['auth_seq_id']
        
        entity = tab['label_asym_id'].unique()
        entity = dict(zip(entity, range(1, len(entity) + 1)))
        tab['label_entity_id'] = tab['label_asym_id'].replace(entity)
        
        tab = tab[
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
        struct = structure(tab)
        struct.name = self.name
        struct.fmt  = self.fmt
        return struct
    
    
    def to_cif(self, path: 'str') -> 'None':
        folder = os.path.split(path)[0]
        if folder:
            os.makedirs(folder, exist_ok=True)
        
        file = open(path, 'w')
        
        struct = self.as_cif()
        
        title = 'data_{}\n'.format(struct.name.upper())
        file.write(title)
        header  = '# \nloop_\n'
        for col in struct.columns:
            header += '_atom_site.{}\n'.format(col)
        file.write(header)
        struct.to_string(file, header=False, index=False)
        file.write('\n# \n')
        
        file.close()
    
    
    def to_pdb(self, path:'str') -> 'None':
        def one_letter_chain_renaming():
            chain_lbl = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
            struct_chain_lbl = tab['auth_asym_id'].astype(str).unique()
            free_chain_lbl  = [c for c in chain_lbl if c not in struct_chain_lbl]
            
            rnm = {}
            cnt = 0
            for lbl in struct_chain_lbl:
                if len(lbl) > 1:
                    try:
                        rnm[lbl] = free_chain_lbl[cnt]
                        cnt += 1
                    except:
                        raise Exception(
                            'free_chain_lbl is over'
                        )
            if cnt:
                tab['auth_asym_id'].replace(rnm, inplace=True)
                for k, v in sorted(rnm.items(), key=lambda x:x[1]):
                    rmk = _REMARK.format(k, v)
                    rmk += ' ' * (80 - len(rmk)) + '\n'
                    msg.append(rmk)
            return msg
        
        folder = os.path.split(path)[0]
        if folder:
            os.makedirs(folder, exist_ok=True)
        
        file = open(path, 'w')
        
        tab = pd.DataFrame(self)
        msg = []
        if self.fmt != 'PDB':
            msg = one_letter_chain_renaming()
        
        text = ''.join(msg)
        for model_num, tt in tab.groupby('pdbx_PDB_model_num', sort=False):
            text += _MODEL.format(model_num)
            tt['id'] = range(1, len(tt) + 1)
            chain_count = 0
            for asym_id, ttt in tt.groupby('auth_asym_id', sort=False):
                ttt['id'] = (ttt['id'] + chain_count) % 1_000_000
                for item in ttt.iloc:
                    text += _ATOM.format(**item)
                item['id'] += 1
                text  += _TER.format(**item)
                chain_count += 1
            text += _ENDMDL
        file.write(text)
    
    def get_mask_by_res(struct:'structure', res:str) -> 'np.ndarray':
        jsls = read_res(res)
        
        if jsls:
            main_mask = np.array([False] * len(struct))
            
            for js in read_res(res):
                cur_mask = np.array([True] * len(struct))
                
                if js['#']:
                    cur_mask &= struct['pdbx_PDB_model_num'].eq(
                        js['#']
                    ).values
                
                if js['/']:
                    cur_mask &= struct['auth_asym_id'].eq(
                        js['/']
                    ).values
                
                if js[':']:
                    if js[':']['B']:
                        cur_mask &= struct['auth_comp_id'].eq(
                            js[':']['B']
                        ).values
                    
                    if js[':']['2']:
                        cur_mask &= struct['auth_seq_id'].between(
                            js[':']['1'], 
                            js[':']['2']
                        ).values
                    
                    else:
                        cur_mask &= struct['auth_seq_id'].eq(
                            js[':']['1']
                        ).values
                
                main_mask |= cur_mask
        else:
            main_mask = np.array([True] * len(struct))
        
        return main_mask
    
    
    def get_substruct_by_res(self, res:'str', 
                             neg:'bool' = False) -> 'structure':
        
        mask = self.get_mask_by_res(res)
        if neg:
            mask ^= True
        struct = structure(self[mask])
        struct.name = self.name
        struct.fmt  = self.fmt
        
        return struct
    
    def describe_residues(self, resrepr:'dict') -> 'list':
        coord = ['Cartn_x', 'Cartn_y', 'Cartn_z']
        
        table = pd.DataFrame(self)
        table.set_index('auth_atom_id', inplace=True)
        table['pdbx_PDB_ins_code'].replace('?', '', inplace=True)
        
        desc = []
        for dssr, res in table.groupby(self.dssr_columns, sort=False):
            flg = False
            res_type = dssr[2]
            if res_type not in resrepr:
                continue
            
            currepr = resrepr[res_type]
            m = []
            for a in currepr:
                try:
                    v = res.loc[a.split(), coord].values
                    m.append(v.mean(axis=0))
                except:
                    flg = True
            if flg:
                continue
            
            m = np.array(m)
            desc.append((dssr, m))
        return desc


def read_cif(path: 'str', name: 'str' = None) -> structure:
    al_col_pairs = (
        ('auth_asym_id', 'label_asym_id'),
        ('auth_seq_id',  'label_seq_id'),
        ('auth_comp_id', 'label_comp_id'),
        ('auth_atom_id', 'label_atom_id'),
    )
    
    if not name:
        name = os.path.splitext(path)[0]
        name = name.split(os.sep)[-1]
    
    with open(path, 'r') as file:
        text = file.read()
    
    start = text.find('_atom_site.')
    if start == -1:
        raise Exception(
            'File {} does not contain _atom_site table.'.format(path)
        )
    end = text.find('#', start) - 1
    tab = text[start:end].split('\n')
    
    columns = []
    for i, line in enumerate(tab):
        if line.startswith('_'):
            columns.append(line.split('.', 1)[1].strip())
        else:
            break
    
    items  = map(str.split, tab[i:])
    struct = structure(items, columns=columns, name=name, fmt='CIF')
    struct = struct.apply(pd.to_numeric)
    
    for a, l in al_col_pairs:
        if a not in struct.columns:
            struct[a] = struct[l]
    
    l = lambda x: x[1:-1] if x.startswith('"') or x.startswith("'") else x
    struct['label_atom_id'] = list(map(l, struct['label_atom_id']))
    struct['auth_atom_id']  = list(map(l, struct['auth_atom_id']))
    
    return struct



def read_pdb(path: 'str', name: 'str' = '') -> 'structure':
    if not name:
        name = os.path.splitext(path)[0]
        name = name.split(os.sep)[-1]
    
    columns = (
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
        
        'pdbx_PDB_model_num',   # extra column for the atomic coordinate table
    )
    
    rec_names = {'ATOM  ', 'HETATM'}
    cur_model = 1
    items     = []
    
    file = open(path, 'r')
    for line in file:
        rec = line[0:6]
        if rec in rec_names:
            item = [
                line[0 : 6].strip(), # group_PDB
                line[6 :11].strip(), # id
                line[12:16].strip(), # auth_atom_id
                line[16:17].strip(), # label_alt_id
                line[17:20].strip(), # auth_comp_id
                line[20:22].strip(), # auth_asym_id
                line[22:26].strip(), # auth_seq_id
                line[26:27].strip(), # pdbx_PDB_ins_code
                line[30:38].strip(), # Cartn_x
                line[38:46].strip(), # Cartn_y
                line[46:54].strip(), # Cartn_z
                line[54:60].strip(), # occupancy
                line[60:66].strip(), # B_iso_or_equiv
                line[76:78].strip(), # type_symbol
                line[78:80].strip(), # pdbx_formal_charge
                cur_model            # pdbx_PDB_model_num
            ]
            items.append(item)
        elif rec == 'MODEL ':
            cur_model = int(line.split()[1])
    file.close()
    
    struct = structure(items, columns=columns, name=name, fmt='PDB')
    struct = struct.apply(pd.to_numeric)
    struct.fillna('', inplace=True)
    
    l = lambda x: x[1:-1] if x.startswith('"') or x.startswith("'") else x
    struct['auth_atom_id']  = list(map(l, struct['auth_atom_id']))
    
    return struct


def read_res(res: str) -> list:
    jsls = []
    subres = res.split()
    for sr in subres:
        js = {
            '#': None, # Model
            '/': None, # Chain
            ':': None  # Residues
        }
        
        # split subres by # / : keys of js
        seq = iter(sr)
        c   = next(seq, False)
        while c:
            if c in js:
                js[c] = ''
                cc = next(seq, False)
                while cc and cc not in js:
                    js[c] += cc
                    cc = next(seq, False)
                else:
                    c = cc
            else:
                c = next(seq, False)
        
        # Transfer string model No to integer
        if js['#']:
            js['#'] = int(js['#'])
        elif js['#'] == None:
            # By default model No = 1
            js['#'] = 1
        
        # Correct / value
        if js['/'] == '':
            js['/'] = None
        
        # Parse : key
        if js[':']:
            elem = js[':'].split('_')
            btp, num1, num2  = None, None, None
            
            btp = elem[0]
            
            if not elem[1].isdigit():
                num1 = ''
                for i, c in enumerate(elem[1]):
                    num1 += c
                num1 = int(num1)
                btp = elem[1][i:]
            else:
                num1 = int(elem[1])
            try:
                num2 = int(elem[2])
            except:
                pass
            
            if btp == '':
                btp = None
            
            js[':'] = {
                'B': btp,
                '1': num1,
                '2': num2 
            }
        
        jsls.append(js)
    
    return jsls



#######
# PDB save formats
_ATOM   = '{group_PDB:<6}{id:>5} {auth_atom_id:<4}{label_alt_id:1}{auth_comp_id:>3}{auth_asym_id:>2}{auth_seq_id:>4}{pdbx_PDB_ins_code:1}   {Cartn_x:>8.3f}{Cartn_y:>8.3f}{Cartn_z:>8.3f}{occupancy:>6.2f}{B_iso_or_equiv:>6.2f}          {type_symbol:>2}{pdbx_formal_charge:>2}\n'
_TER    = 'TER   {id:>5}      {auth_comp_id:>3}{auth_asym_id:>2}{auth_seq_id:>4}                                                      \n'
_MODEL  = 'MODEL     {:>4}                                                                  \n'
_REMARK = 'REMARK 250 CHAIN RENAMING {} -> {}'
_ENDMDL = 'ENDMDL                                                                          \n'
