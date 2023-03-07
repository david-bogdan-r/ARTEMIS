import os
import pandas as pd
import numpy  as np

from copy   import deepcopy
from string import ascii_letters, digits


pd.to_numeric.__defaults__ = 'ignore', None

COORD = ['Cartn_x', 'Cartn_y', 'Cartn_z']

class Structure:
    def __init__(self, name:'str', atom_site:'pd.DataFrame'):
        self.name      = name
        self.atom_site = atom_site
    
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return '<{} Structure>'.format(self)
    
    def to_cif(self, path:'str'):
        folder, file = os.path.split(path)
        if folder:
            os.makedirs(folder, exist_ok=True)
        
        with open(path, 'w') as file:
            file.write('data_{}\n'.format(self))
            file.write('# \nloop_\n')
            file.write('\n'.join(self.atom_site.add_prefix('_atom_site.')) + '\n')
            self.atom_site.to_string(file, header=False, index=False)
            file.write('\n# \n')
    
    
    def to_pdb(self, path:'str'):
        folder, file = os.path.split(path)
        if folder:
            os.makedirs(folder, exist_ok=True)
        
        atom_site = self.atom_site.copy()
        atom_site.replace('.', '', inplace=True)
        atom_site.replace('?', '', inplace=True)
        
        with open(path, 'w') as file:
            chains = set(self.atom_site['auth_asym_id'].astype(str))
            rename = {c: None for c in chains if len(c) > 1}
            if rename:
                names = (name for name in ascii_letters + digits)
                for c in rename:
                    name = next(names, None)
                    while name in chains:
                        name = next(names, None)
                    if name:
                        rename[c] = name
                    else:
                        raise ValueError('No free names for the chain')
                    
                    file.write(PDB_FORMAT['REMARK'].format(c, rename[c]))
                atom_site = atom_site.replace({'auth_asym_id':rename})
            
            chain_delta = 0
            model_delta = 0
            
            model = atom_site.iloc[0]['pdbx_PDB_model_num']
            chain = atom_site.iloc[0]['auth_asym_id']
            
            file.write(PDB_FORMAT['MODEL'].format(model))
            item_ = None
            for item in atom_site.to_dict('index').values():
                if item['pdbx_PDB_model_num'] != model:
                    model = item['pdbx_PDB_model_num']
                    file.write(PDB_FORMAT['ENDMDL'])
                    file.write(PDB_FORMAT['MODEL'].format(model))
                    model_delta += 1
                elif item['auth_asym_id'] != chain:
                    chain = item['auth_asym_id']
                    item_['id'] += 1
                    file.write(PDB_FORMAT['TER'].format(**item_))
                    chain_delta += 1
                
                item['id'] = (item['id'] + chain_delta + model_delta) % 100_000
                item_ = item
                file.write(PDB_FORMAT['ATOM'].format(**item))
    
    def get_coord(self):
        return self.atom_site[COORD].values
    
    def set_coord(self, mat:'np.ndarray'):
        self.atom_site.loc[:, COORD] = mat
    
    coord = property(get_coord, set_coord)
    
    def copy(self):
        return deepcopy(self)
    
    def __add__(self, val):
        other = self.copy()
        other.coord +=  val
        return other
    
    def __sub__(self, val):
        other = self.copy()
        other.coord -=  val
        return other
    
    def __mul__(self, val):
        other = self.copy()
        if hasattr(val, 'shape') and val.shape == (3, 3):
            other.coord = np.dot(other.coord, val)
        else:
            other.coord *= val
        return other
    
    @property
    def dssr(self):
        columns = [
            'pdbx_PDB_model_num',
            'auth_asym_id',
            'auth_comp_id',
            'auth_seq_id',
            'pdbx_PDB_ins_code'
        ]
        return self.atom_site[columns].replace('?', '').values
    
    def drop_altloc(self, keep='last'):
        self.atom_site.drop_duplicates(
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
    
    def _read_res(self, res:'str'):
        reads  = []
        for subres in res.split():
            read = {
                '#': None,
                '/': None,
                ':': None
            }
            
            ic = iter(subres)
            c  = next(ic, None)
            while c:
                if c in read:
                    read[c] = ''
                    cc = next(ic, None)
                    while cc and cc not in read:
                        read[c] += cc
                        cc = next(ic, None)
                    else:
                        c = cc
                else:
                    c = next(ic, None)
            
            if read['#']:
                read['#'] = int(read['#'])
            
            if read['/'] == '':
                read['/'] = None
            
            if read[':']:
                part      = read[':'].split('_')
                b, n1, n2 = None, None, None
                
                if len(part) == 1:
                    b = part[0]
                
                elif len(part) == 2:
                    if part[0]:
                        b = part[0]
                    
                    if part[1].startswith('-'):
                        n1 = '-'
                        part[1] = part[1][1:]
                    else:
                        n1 = ''
                    
                    if not part[1].isdigit():
                        for i, c in enumerate(part[1]):
                            if c.isdigit():
                                n1 += c
                            else:
                                break
                        n1 = int(n1)
                        b  = part[1][i:]
                    else:
                        n1 = int(n1 + part[1])
                
                elif len(part) == 3:
                    b = part[0]
                    n1 = int(part[1])
                    n2 = int(part[2])
                
                read[':'] = {'b':b, 'n1':n1, 'n2':n2}
            reads.append(read)
        
        return reads
    
    def _res_mask(self, res:'str'):
        reads = self._read_res(res)
        if reads:
            mask = np.array([False]*len(self.atom_site))
            for read in reads:
                m = np.array([True]*len(self.atom_site))
                if read['#']:
                    m &= (self.atom_site['pdbx_PDB_model_num']
                          .eq(read['#']).values)
                if read['/']:
                    m &= (self.atom_site['auth_asym_id']
                          .eq(read['/']).values)
                if read[':']:
                    if read[':']['b']:
                        m &= (self.atom_site['auth_comp_id']
                              .eq(read[':']['b']).values)
                    if read[':']['n2']:
                        m &= (self.atom_site['auth_seq_id']
                              .between(read[':']['n1'], read[':']['n2']).values)
                    else:
                        m &= (self.atom_site['auth_seq_id']
                              .eq(read[':']['n1']).values)
                mask |= m
        else:
            mask = np.array([True]*len(self.atom_site))
        
        return mask
    
    def get_substruct(self, res:'str', neg:'bool'=False):
        mask = self._res_mask(res)
        if neg:
            mask ^= True
        structure = Structure(self.name, self.atom_site[mask])
        return structure
    
    def resgen(self, res:'str'='', neg:'bool'=False):
        substruct = self.get_substruct(res, neg)
        
        columns = [
            'pdbx_PDB_model_num',
            'auth_asym_id',
            'auth_comp_id',
            'auth_seq_id',
            'pdbx_PDB_ins_code'
        ]
        
        for dssr, atom_site in substruct.atom_site.groupby(columns, sort=False):
            name = self.name + '_' + '.'.join(map(str, dssr)).replace('?', '')
            res  = Structure(name, atom_site)
            yield res
    
    def resapply(self, func):
        results = []
        for res in self.resgen(''):
            results.append(func(res))
        return results


def resrepr(res:'Structure', rep:'dict'={}):
    base = res.atom_site['auth_comp_id'].iloc[0]
    if base not in rep:
        return None
    else:
        atom_site:'pd.DataFrame' = res.atom_site.copy()
        atom_site.set_index('auth_atom_id', inplace=True)
        
        mat = []
        for atom in rep[base]:
            try:
                vec = atom_site.loc[atom.split(), COORD].values
                mat.append(vec.mean(axis=0))
            except:
                return None
        return np.array(mat)


def join(*structure, name:'str'='join',):
    atom_site = pd.concat([s.atom_site for s in structure])
    structure = Structure(name, atom_site)
    return structure


def read_pdb(path:'str', name:'str'=''):
    if not name:
        name = os.path.basename(path).split(os.path.extsep)[0]
    
    index = {'ATOM  ', 'HETATM'}
    model = 1
    items = []
    for line in open(path, 'r'):
        ind = line[:6]
        if ind in index:
            item = [
                line[0 : 6],
                line[6 :11],
                line[12:16],
                line[16:17],
                line[17:20],
                line[20:22],
                line[22:26],
                line[26:27],
                line[30:38],
                line[38:46],
                line[46:54],
                line[54:60],
                line[60:66],
                line[76:78],
                line[78:80]
            ]
            item = list(map(str.strip, item))
            
            auth_atom_id = item[2]
            if auth_atom_id.startswith('"') or auth_atom_id.startswith("'"):
                item[2] = auth_atom_id[1:-1]
            
            item += [model]
            items.append(item)
        elif ind == 'MODEL ':
            model = int(line.split()[1])
    
    atom_site = pd.DataFrame(
        items,
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
    
    atom_site['label_atom_id'] = atom_site['auth_atom_id']
    atom_site['label_comp_id'] = atom_site['auth_comp_id']
    atom_site['label_asym_id'] = atom_site['auth_asym_id']
    atom_site['label_seq_id']  = atom_site['auth_seq_id']
    
    entity_id = atom_site['label_asym_id'].unique()
    entity_id = dict(zip(entity_id, range(1, len(entity_id) + 1)))
    atom_site['label_entity_id'] = atom_site['label_asym_id'].replace(entity_id)
    
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
    
    structure = Structure(name, atom_site)
    
    return structure


def read_cif(path:'str', name:'str'=''):
    if not name:
        name = os.path.basename(path).split(os.path.extsep)[0]
    
    with open(path, 'r') as file:
        data = file.read()
    
    start = data.find('_atom_site.')
    end   = data.find('#', start) - 1
    data  = data[start:end].split('\n')
    
    columns = []
    for i, line in enumerate(data):
        if line.startswith('_'):
            columns.append(line.split('.')[1].strip())
        else:
            break
    items = map(str.split, data[i:])
    atom_site = pd.DataFrame(items, columns=columns)
    
    quotes = atom_site['label_atom_id'].str.startswith(('"', '"'))
    if quotes.any():
        atom_site.loc[quotes, 'label_atom_id'] = \
            atom_site[quotes]['label_atom_id'].str[1:-1]
    if 'auth_atom_id' in atom_site:
        quotes = atom_site['auth_atom_id'].str.startswith(('"', '"'))
        if quotes.any():
            atom_site.loc[quotes, 'auth_atom_id'] = \
                atom_site[quotes]['auth_atom_id'].str[1:-1]
    
    atom_site = atom_site.apply(pd.to_numeric)
    structure = Structure(name, atom_site)
    
    return structure


PDB_FORMAT = {
    'ATOM'  : '{group_PDB:<6}{id:>5} {auth_atom_id:<4}{label_alt_id:1}{auth_comp_id:>3}{auth_asym_id:>2}{auth_seq_id:>4}{pdbx_PDB_ins_code:1}   {Cartn_x:>8.3f}{Cartn_y:>8.3f}{Cartn_z:>8.3f}{occupancy:>6.2f}{B_iso_or_equiv:>6.2f}          {type_symbol:>2}{pdbx_formal_charge:>2}\n',
    'TER'   : 'TER   {id:>5}      {auth_comp_id:>3}{auth_asym_id:>2}{auth_seq_id:>4}                                                      \n',
    'MODEL' : 'MODEL     {:>4}                                                                  \n',
    'REMARK': 'REMARK 250 CHAIN RENAMING {} -> {}\n',
    'ENDMDL': 'ENDMDL                                                                          \n',
}
