from ase.build import molecule
from ase import Atoms
from ase.calculators.dmol import DMol3
import numpy as np
from pathlib import Path

class DMol3KudoCop():
    """DMol3を用いて第一原理計算するクラス
    """
    def __init__():
        pass

    def dmol3_calc(self, calc_directory='dmol_calc', **kwargs):
        positions = self[['x','y','z']]
        symbols = self['type'].map(self.atom_type_to_symbol)
        cell = self.cell
        self.atoms_calc = Atoms(
            positions=positions,
            symbols=symbols,
            cell=cell)
        calc = DMol3(**kwargs)
        calc.directory = calc_directory
        self.atoms_calc.calc = calc
    
