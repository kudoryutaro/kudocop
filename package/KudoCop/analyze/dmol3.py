from ase.build import molecule
from ase import Atoms
from ase.calculators.dmol import DMol3
import numpy as np
from pathlib import Path
from tqdm import trange, tqdm

class DMol3KudoCop():
    """DMol3を用いて第一原理計算するクラス
    """
    def __init__():
        pass

    def dmol3_calc(self, calc_directory='dmol_calc', **kwargs):
        """dmol3で第一原理計算をする関数
        Parameters
        ----------
            calc_directory : Path
                このフォルダ内にdmol3の計算した結果を保存する

        Note
        ----
            dmol3_calcを用いて第一原理計算をしたあとにtotal_energyやforcesなどを取り出すことが出来る
            dmol3_calcを行った後にatoms_calcというプロパティができるのでそこからenergyなどを取り出す

            Examples
            --------
                >>>sdat.atoms_calc.get_total_energy() # 単位:eV
                -1537.0315543

                >>>sdat.atoms_calc.get_forces() # 単位:eV / Å
                [[ 0.00000000e+00  0.00000000e+00  1.62750973e+00]
                [ 1.34485812e+00  2.32936367e+00 -5.42503242e-01]
                [-2.68971624e+00  4.57251760e-09 -5.42503242e-01]
                [ 1.34485811e+00 -2.32936368e+00 -5.42503242e-01]]
        """
        calc_directory = Path(calc_directory)
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
    

class DMol3KudoCopForSDats():
    """DMol3を用いて第一原理計算するクラス
    """
    def __init__():
        pass

    def dmol3_calc(self, calc_directory='dmol_calc', **kwargs):
        """dmol3で第一原理計算をする関数
        Parameters
        ----------
            calc_directory : Path
                このフォルダ内にdmol3の計算した結果を保存する
        Note
        ----
            dmol3_calcを用いて第一原理計算をしたあとにtotal_energyやforcesなどを取り出すことが出来る
            dmol3_calcを行った後にatoms_calcというリストのプロパティができる
            sdat.atoms_calc[step_idx].メソッド() で個別のenergyなどを取り出す
            または、sdat.get_total_energies()などで、リスト形式ですべて取り出す
            Examples
            --------
                step_idx = 0
                >>>sdat.atoms_calc[step_idx].get_total_energy() # 単位:eV
                -1537.0315543

                >>>sdat.atoms_calc[step_idx].get_forces() # 単位:eV / Å
                [[ 0.00000000e+00  0.00000000e+00  1.62750973e+00]
                [ 1.34485812e+00  2.32936367e+00 -5.42503242e-01]
                [-2.68971624e+00  4.57251760e-09 -5.42503242e-01]
                [ 1.34485811e+00 -2.32936368e+00 -5.42503242e-01]]
        """
        calc_directory = Path(calc_directory)
        self.atoms_calc = [None for _ in range(len(self.step_nums))]

        for step_idx in range(len(self.step_nums)):
            positions = self.atoms[step_idx][['x','y','z']]
            symbols = self.atoms[step_idx]['type'].map(self.atom_type_to_symbol)
            cell = self.cell
            self.atoms_calc[step_idx] = Atoms(
                positions=positions,
                symbols=symbols,
                cell=cell)
            calc = DMol3(**kwargs)
            calc.directory = calc_directory / f'calc_step_idx_{step_idx}'
            self.atoms_calc[step_idx].calc = calc
    
    def get_dmol3_total_energies(self):
        """計算したすべてのフレームごとのエネルギーをリスト形式で返す関数
        Returns:
            total_energies : list
                total_energies[step_idx]でstep_idxのエネルギーを取り出すことができる
        """
        assert hasattr(self, 'atoms_calc'), 'dmol3_calcで先に第一原理計算してください'

        total_energies = [None] * len(self.step_nums)
        for step_idx in trange(len(self.step_nums)):
            total_energies[step_idx] = self.atoms_calc[step_idx].get_total_energy()
        return total_energies
    
    
    def get_dmol3_forces(self):
        """計算したすべてのフレームごとのforceをリスト形式で返す関数
        Returns:
            total_energies : list
                total_energies[step_idx]でstep_idxのforceを取り出すことができる
        """
        assert hasattr(self, 'atoms_calc'), 'dmol3_calcで先に第一原理計算してください'

        forces = [None] * len(self.step_nums)
        for step_idx in trange(len(self.step_nums)):
            forces[step_idx] = self.atoms_calc[step_idx].get_forces()
        return forces
    

