import numpy as np
from pathlib import Path
from tqdm import trange, tqdm
import subprocess
from pathlib import Path
import os
import time

try:
    from ase.build import molecule
    from ase import Atoms
    from ase.calculators.dmol import DMol3
except:
    pass

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

    def dmol3_md(self, calc_label='dmol3_md', calc_directory='dmol3_md', np=1,
                ensemble='NVE', temperature=300.0, time_step=1.0, number_of_steps=1000, exist_ok=False,
                max_memory=2048, print_outmol=True):
        """DMol3を用いて第一原理分子動力学を行う.

        Parameters
        ----------
        calc_label : str
            計算する名前 計算結果は calc_directory/calc_label.outmol になる
        calc_directory : str or Path
            計算するディレクトリ
        np : int
            並列計算数
        ensemble : str ['NVE', 'NVT']
            アンサンブルの種類
        temperature : float
            温度, 単位は[K]
        time_step : float
            タイムステップ, 単位は[fs]
        number_of_steps : int
            ステップ数
        exist_ok : bool
            exist_ok = Trueとすると、フォルダがあっても上書きする
            exist_ok = Falseとすると、フォルダが有るときは上書きしない
            デフォルトはFalse
        max_memory : int
            使用する最大のメモリ
        print_outmol : bool
            実行中にoutmolファイルを標準出力するかどうか
        """    

        assert ensemble in ['NVE', 'NVT'], f"ensemble : {ensemble}には対応していません. ['NVE', 'NVT']のみ対応"
        
        calc_directory = Path(calc_directory)
        os.makedirs(calc_directory, exist_ok=exist_ok)
        dmol3_input_path = calc_directory / f'{calc_label}.input'
        dmol3_input_lines = f"""

# Task parameters
Calculate                     molecular_dynamics
Write_HIS_File                on
Write_ARC_File                on
MD_Velocity                   random
MD_Time_Step                  {time_step:.4f}
MD_Simann_panel 
{int(number_of_steps)}      MD_{ensemble}     {temperature:.4f}
# 
Symmetry                      off
Max_memory                    {int(max_memory)}
File_usage                    smart
Scf_density_convergence       1.000000e-06
Scf_charge_mixing             2.000000e-01
Scf_diis                      6 pulay
Scf_iterations                50

# Electronic parameters
Spin_polarization             restricted
Charge                        0
Forces                        on
Basis                         dnp
Pseudopotential               none
Functional                    pwc
Harris                        off
Aux_density                   hexadecapole
Integration_grid              fine
Occupation                    fermi
Cutoff_Global                 3.4000 angstrom

# Calculated properties

        """
        with open(dmol3_input_path, 'w') as f:
            f.write(dmol3_input_lines)
        self.export_car(calc_directory / f'{calc_label}.car')
        cmd = f'RunDMol3.sh {calc_label} -np {np}'
        dmol_md_process = subprocess.Popen(cmd, cwd=calc_directory, shell=True)
        if print_outmol:
            tail_process = subprocess.Popen(f'tail -F {calc_label}.outmol', cwd=calc_directory, shell=True)
            while dmol_md_process.poll() is None:
                time.sleep(1)
            tail_process.kill()


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
    

