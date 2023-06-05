import numpy as np
import pandas as pd
from tqdm import trange
from pathlib import Path
import subprocess
import os
import time

class Laich():
    """laichを用いて分子動力学計算をするクラス
    """
    def __init__():
        pass

    def laich_md(self, calc_directory='laich_calc', para_file_path=None, laich_cmd='laich',
                time_step=0.25, total_step=10000, file_step=1000, save_restart_step=10000, mpi_grid_x=1, mpi_grid_y=1, mpi_grid_z=1, omp_grid_x=1, omp_grid_y=1, omp_grid_z=1, cut_off=10.0, margin=1.0, 
                ghost_factor=20.0, show_mask=1, read_velocity=0, thermo='Langevin',aim_temp=300.0, 
                init_temp=300.0, final_temp=300.0, thermo_freq=0.005, 
                sel=50, weight_path='./script_model.pth', ngpus=1, nnp_model_path='script_model.pth', 
                exist_ok=False):
        """laichで分子動力学を実行する関数
            calc_directory : str or Path
                laichで構造最適化をする時に使用するディレクトリ
            para_file_path : str or Path
                laichで構造最適化する時に使用するReaxFFのパラメータ
            laich_cmd : str
                laichを実行するコマンド
            exist_ok : bool
                exist_ok = Trueのとき、calc_directoryが存在していてもエラーは出ずに上書きする
                exist_ok = Falseのとき、calc_directoryが存在しているときはエラーが出る
            time_step : float
            total_step : int
            file_step : int
            save_restart_step : int
            mpi_grid_x : int
            mpi_grid_y : int
            mpi_grid_z : int
            cut_off : float
            margin : float
            ghpst_factor : float
            show_mask : float
            read_velocity : int
            thremo : str
            aim_temp : float
            init_temp : float
            final_tmep : float
            thermo_freq : float
                Laichのドキュメントを参照してください
        """
        calc_directory = Path(calc_directory) 
        if para_file_path is not None:
            para_file_path = Path(para_file_path)
        input_file_path = calc_directory / 'input.rd'
        config_file_path = calc_directory / 'config.rd'
        laich_md_config = f"""Mode MD
ForceField Reaxff
XYZFile input.rd
ParaFile para.rd
TimeStep {time_step}
TotalStep {total_step}
ObserveStep 1
FileStep {file_step}
BondStep {file_step}
SaveRestartStep {save_restart_step}
SEL {sel}
NPUGS {ngpus}
NNPModelPath {nnp_model_path}
WEIGHTPATH {weight_path}
MPIGridX {mpi_grid_x}
MPIGridY {mpi_grid_y}
MPIGridZ {mpi_grid_z}
OMPGridX {omp_grid_x}
OMPGridY {omp_grid_y}
OMPGridZ {omp_grid_z}
CUTOFF {cut_off}
MARGIN {margin}
ShowMask {show_mask}
ReadVelocity {read_velocity}
Thermo  {thermo}
AimTemp {aim_temp}
InitTemp {init_temp}
FinalTemp {final_temp}
ThermoFreq  {thermo_freq}
GhostFactor {ghost_factor}
        """
        os.makedirs(calc_directory, exist_ok=exist_ok)
        self.export_input(input_file_path)
        with open(config_file_path, 'w') as f:
            f.writelines(laich_md_config)
            
        if para_file_path is not None:
            subprocess.run(f'cp {para_file_path} {calc_directory / "para.rd"}', shell=True)
        
        np = mpi_grid_x * mpi_grid_y * mpi_grid_z

        cmd = f"mpiexec.hydra -np {np} {laich_cmd} < /dev/null >& out"
        laich_md_process = subprocess.Popen(cmd, cwd=calc_directory, shell=True)

        laich_md_process.wait()

        dumppos_file_paths = list(calc_directory.glob('./dump.pos.*'))
        dumppos_file_paths.sort(reverse=True)
        dumpbond_file_paths = list(calc_directory.glob('./dump.bond.*'))
        dumpbond_file_paths.sort(reverse=True)
        if len(dumppos_file_paths) == 0:
            raise RuntimeError('dumpposが生成されていません')
        if len(dumpbond_file_paths) == 0:
            raise RuntimeError('dumpbondが生成されていません')
        optimized_dumppos_file_path = dumppos_file_paths[0]
        self.import_dumppos(optimized_dumppos_file_path)
        optimized_dumpbond_file_path = dumpbond_file_paths[0]
        self.import_dumpbond(optimized_dumpbond_file_path)

    def laich_opt(self, calc_directory='laich_calc', para_file_path=None, laich_cmd='laich',
                time_step=0.25, total_step=10000, file_step=1000, save_restart_step=10000, mpi_grid_x=1, mpi_grid_y=1, mpi_grid_z=1, cut_off=10.0, margin=1.0, 
                sel=50, weight_path='./script_model.pth', ngpus=1, 
                ghost_factor=20.0, del_r=0.0001, max_r=0.1, exist_ok=False):
        """laichで構造最適化する関数
        Parameters
        ----------
            calc_directory : str or Path
                laichで構造最適化をする時に使用するディレクトリ
            para_file_path : str or Path
                laichで構造最適化する時に使用するReaxFFのパラメータ
            laich_cmd : str
                laichを実行するコマンド
            exist_ok : bool
                exist_ok = Trueのとき、calc_directoryが存在していてもエラーは出ずに上書きする
                exist_ok = Falseのとき、calc_directoryが存在しているときはエラーが出る
            time_step : float
            total_step : int
            file_step : int
            save_restart_step : int
            mpi_grid_x : int
            mpi_grid_y : int
            mpi_grid_z : int
            cut_off : float
            margin : float
            ghpst_factor : float
            del_r : float
            max_r float
                Laichのドキュメントを参照してください
        """
        calc_directory = Path(calc_directory) 
        if para_file_path is not None:
            para_file_path = Path(para_file_path)
        input_file_path = calc_directory / 'input.rd'
        config_file_path = calc_directory / 'config.rd'
        lacih_opt_config = f"""Mode OPT
ForceField Reaxff
XYZFile input.rd
ParaFile para.rd
TimeStep {time_step}
TotalStep {total_step}
ObserveStep 1
FileStep {file_step}
BondStep {file_step}
SaveRestartStep {save_restart_step}
MPIGridX {mpi_grid_x}
MPIGridY {mpi_grid_y}
MPIGridZ {mpi_grid_z}
SEL {sel}
NPUGS {ngpus}
WEIGHTPATH {weight_path}
CUTOFF {cut_off}
MARGIN {margin}
GhostFactor {ghost_factor}
DelR {del_r}
MaxR {max_r}
""" 
        os.makedirs(calc_directory, exist_ok=exist_ok)
        self.export_input(input_file_path)
        with open(config_file_path, 'w') as f:
            f.writelines(lacih_opt_config)
        
        if para_file_path is not None:
            try:
                subprocess.run(f'cp {para_file_path} {calc_directory / "para.rd"}', shell=True)
            except:
                pass
        
        np = mpi_grid_x * mpi_grid_y * mpi_grid_z

        cmd = f"mpiexec.hydra -np {np} {laich_cmd} < /dev/null >& out"
        laich_opt_process = subprocess.Popen(cmd, cwd=calc_directory, shell=True)

        laich_opt_process.wait()

        dumppos_file_paths = list(calc_directory.glob('./dump.pos.*'))
        dumppos_file_paths.sort(reverse=True)
        dumpbond_file_paths = list(calc_directory.glob('./dump.bond.*'))
        dumpbond_file_paths.sort(reverse=True)
        if len(dumppos_file_paths) == 0:
            raise RuntimeError('dumpposが生成されていません')
        if len(dumpbond_file_paths) == 0:
            raise RuntimeError('dumpbondが生成されていません')
        optimized_dumppos_file_path = dumppos_file_paths[0]
        self.import_dumppos(optimized_dumppos_file_path)
        optimized_dumpbond_file_path = dumpbond_file_paths[0]
        self.import_dumpbond(optimized_dumpbond_file_path)