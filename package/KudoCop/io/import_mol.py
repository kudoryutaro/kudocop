import numpy as np
import pandas as pd
import sys
from pathlib import Path

try:
    from ase import build
except:
    pass

class ImportMol():
    def __init__(self):
        pass

    def import_mol(self, molecular_fomula):
        '''分子式から原子配置を読み込む
        Parameters
        ----------
            molecular_fomula : str
                読み込む分子式
        '''
        ase_atoms = build.molecule(molecular_fomula)
        self.atoms = pd.DataFrame(ase_atoms.positions, columns=['x', 'y', 'z'])
        self.atoms['type'] = ase_atoms.get_chemical_symbols()
        self.atoms['type'] = self.atoms['type'].map(self.atom_symbol_to_type)
        self.atoms['x'] += abs(self.atoms['x'].min()) + 0.1
        self.atoms['y'] += abs(self.atoms['y'].min()) + 0.1
        self.atoms['z'] += abs(self.atoms['z'].min()) + 0.1
        