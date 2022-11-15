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
                ensemble='NVE', temperature=300.0, time_step=1.0, number_of_steps=1000, exist_ok=False, max_memory=2048, print_outmol=True, scf_density_convergence=1.000000e-05, 
                integration_grid='medium', basis='dnp', cutoff_Global=3.2000, scf_iterations=50, run=False, save_rundmol3_sh=False, save_qsub_script=False):
        """DMol3を用いて第一原理分子動力学を行う.

        Parameters
        ----------
        calc_label : str
            計算する名前 計算結果は calc_directory/calc_label.outmol になる
        calc_directory : str or Path
            計算するディレクトリ
        np : int
            並列計算数
        run : bool
            run=Trueのときは実際にdmol3で第一原理計算MDを実行する
            run=Falseのときは設定ファイルを書き込むだけ
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
        scf_density_convergence : float
        integration_grid : str ['coarse', 'medium', 'fine']
        basis : str  ['dn', 'dnd', 'dnp']
        cutoff_Global : float
        scf_iterations : int
            DMol3のマニュアルを参照してください
        """    

        assert ensemble in ['NVE', 'NVT'], f"ensemble : {ensemble}には対応していません. ['NVE', 'NVT']のみ対応"

        calc_directory = Path(calc_directory)
        os.makedirs(calc_directory, exist_ok=exist_ok)
        dmol3_input_path = calc_directory / f'{calc_label}.input'

        if ensemble == 'NVE':
            MD_Simann_panel = f'{int(number_of_steps)}      MD_NVE     {temperature:.4f}'
        elif ensemble == 'NVT':
            relaxation_time = 10.0
            chain_length = 2
            MD_Simann_panel = f'   {int(number_of_steps)}      NVT_MGGMT     {temperature:.4f}   {relaxation_time}   {chain_length}'

        
        dmol3_input_lines = f"""

# Task parameters
Calculate                     molecular_dynamics
Write_HIS_File                on
Write_ARC_File                on
MD_Velocity                   random
MD_Time_Step                  {time_step:.4f}
MD_Simann_panel 
{MD_Simann_panel}
# 
Symmetry                      off
Max_memory                    {int(max_memory)}
File_usage                    smart
Scf_density_convergence       {scf_density_convergence:.6e}
Scf_charge_mixing             2.000000e-01
Scf_diis                      6 pulay
Scf_iterations                {int(scf_iterations)}

# Electronic parameters
Spin_polarization             restricted
Charge                        0
Forces                        on
Basis                         {basis}
Pseudopotential               none
Functional                    pwc
Harris                        off
Aux_density                   hexadecapole
Integration_grid              {integration_grid}
Occupation                    fermi
Cutoff_Global                 {cutoff_Global:.4f} angstrom

# Calculated properties

        """
        with open(dmol3_input_path, 'w') as f:
            f.write(dmol3_input_lines)
        self.export_car(calc_directory / f'{calc_label}.car')
        with open(calc_directory / 'RunDMol3.sh', 'w') as f:
            f.write(RunDMol3_sh)
        
        with open(calc_directory / 'qsub_script.sh', 'w') as f:
            qsub_script = get_qsub_script(calc_label)
            f.write(qsub_script)

        if not run:
            return
            
        cmd = f'RunDMol3.sh {calc_label} -np {np}'
        dmol_md_process = subprocess.Popen(cmd, cwd=calc_directory, shell=True)
        time.sleep(5)
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
                cell=cell,
                pbc=[1, 1, 1])
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
    

RunDMol3_sh = """#!/bin/sh
# ---------------------------------------------------------------------------
#
# Script for stand-alone DMol3 execution.
#
# This program and all subroutines, data, and files used by it
# are protected by Copyright and hence may not be used, copied,
# modified, transmitted, inspected, or executed by any means including
# the use of electronic data processing equipment, xerography, or
# any other methods without the express written permission of the
# copyright holder.
#
# Copyright (c) 2015, Dassault Systemes, All Rights Reserved
#
# ***************************************************************************
MS_INSTALL_ROOT=/work/app/MaterialsStudio2021HF1/MaterialsStudio21.1
export MS_INSTALL_ROOT
server=DMol3
$MS_INSTALL_ROOT/share/bin/runMSserver.sh $server "$@"
if [ $? != 6 ]; then
    exit $?
fi

cat $MS_INSTALL_ROOT/etc/DMol3/bin/RunDMol3.Readme"""

def get_qsub_script(calc_label):
    qsub_script = f"""#!/bin/sh
#PBS -l select=1
#PBS -l dmol3=1
#PBS -q C_002
#PBS -N {calc_label}

DIRNAME=`basename $PBS_O_WORKDIR`
WORKDIR=/work/$USER/$PBS_JOBID
mkdir -p $WORKDIR
cp -raf  $PBS_O_WORKDIR $WORKDIR
cd $WORKDIR/$DIRNAME

./RunDMol3.sh -np 18 {calc_label}

cd; if cp -raf $WORKDIR/$DIRNAME $PBS_O_WORKDIR/.. ; then rm -rf $WORKDIR; fi
    """
    return qsub_script
