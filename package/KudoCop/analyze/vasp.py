import numpy as np
import pathlib
from tqdm import trange, tqdm
import subprocess
from typing import List
from pathlib import Path
import os
import time


try:
    from ase.build import molecule
    from ase import Atoms
    from ase.calculators.dmol import DMol3
except:
    pass

class Vasp():
    """vaspを用いて第一原理計算するクラス
    """
    def __init__():
        pass

    def vasp(
            self,
            calc_directory, 
            poscar_ofn: str="POSCAR", 
            poscar_comment: str="", 
            poscar_scaling_factor: float=1.0,
            incar_ofn: str="INCAR",
            incar_config: dict=None,  
            potcar_ofn: str="POTCAR",
            potcar_root: str=None,  
            kpoints_ofn: str="KPOINTS",
            kpoints_comment: str="",
            kpoints_kx: int=1,
            kpoints_ky: int=1,
            kpoints_kz: int=1,
            iconst_ofn: str="ICONST",
            iconst_config: List[str]=None,
            mpi_command: str="aprun",
            vasp_command: str="vasp_std",
            print_vasp: bool=True,
            exist_ok: bool=False
        ):
        """vaspを実行する.
        Parameters
        ----------
            calc_directory: str
                vaspで計算を実行するディレクトリ
            poscar_ofn: str
                POSCARのファイル名, 基本的に変える必要はない
            poscar_comment: str
                POSCARの1行目に書かれるコメント
            poscar_scaling_factor: float
                vaspのドキュメンをを参照してください.基本的に変える必要はない
            incar_ofn: str
                INCARのファイル名, 基本的に変える必要はない
            potcar_root: str
                元のPOTCARがあるフォルダ
            kpoints_ofn: str
                KPOINTSのファイル名, 基本的に変える必要はない
            kpoints_comment: str
                KPOINTSの1行目に書かれるコメント
            kpoints_kx: int
            kpoints_ky: int
            kpoints_kz: int
                x, y, z方向のK点     
            iconst_ofn: str
                ICONSTのファイル名, 基本的に変える必要はない
            iconst_config: List[str]
                ICONSTに書かれるconfig
                VASPでNPT計算をする時で、セルの角度を固定したい時は以下のように指定する
                iconst_config = ['LA 1 2 0',
                            'LA 1 3 0',
                            'LA 2 3 0']
            mpi_command: str
                実行するmpiのコマンド、例えばmpirunなど
            vasp_command: str
                実行するvaspのコマンド
            print_vasp: bool
                実行中のvaspの標準出力をpythonからも表示するかどうか
            exist_ok: bool
                calc_directoryが存在する時、exist_ok=Trueのときは上書きしてvaspを実行する
                calc_directoryが存在する時、exist_ok=Falseのときは上書きしない
        
        """
        assert potcar_root is not None
        assert incar_config is not None
        if "NCORE" not in incar_config:
            incar_config["NCORE"] = 1
        if "NPAR" not in incar_config:
            incar_config["NPAR"] = 1
        if "KPAR" not in incar_config:
            incar_config["KPAR"] = 1
        calc_directory = pathlib.Path(calc_directory)
        os.makedirs(calc_directory, exist_ok=exist_ok)
        #poscar_path = calc_directory / poscar_ofn
        incar_path = calc_directory / incar_ofn
        potcar_path = calc_directory / potcar_ofn
        kpoints_path = calc_directory / kpoints_ofn
        iconst_path = calc_directory / iconst_ofn
        num_process = incar_config["NCORE"] * incar_config["NPAR"] * incar_config["KPAR"]
        #self.export_vasp_poscar(poscar_path, poscar_comment, poscar_scaling_factor)
        self.export_vasp_incar(incar_path, incar_config)
        self.export_vasp_potcar(potcar_path, potcar_root)
        self.export_vasp_kpoints(kpoints_path, kpoints_comment, 
                                 kpoints_kx, kpoints_ky, kpoints_kz)
        if iconst_config is not None:
            self.export_vasp_iconst(iconst_path, iconst_config)

        # For kbox
        #cmd = f'mpiexec -np {num_process} {vasp_command} > stdout'
        # For MASAMUNE-IMR
        if num_process <= 36:
            cmd = f'{mpi_command} -n {num_process} -N {num_process} -j 1 {vasp_command} > stdout'
        else:
            cmd = f'{mpi_command} -n {num_process} -N 36 -j 1 {vasp_command} > stdout'
        dmol_md_process = subprocess.Popen(cmd, cwd=calc_directory, shell=True)
        time.sleep(5)
        if print_vasp:
            tail_process = subprocess.Popen(f'tail -F stdout', cwd=calc_directory, shell=True)
            while dmol_md_process.poll() is None:
                time.sleep(1)
            tail_process.kill()
