import numpy as np
import pandas as pd
import sys
from pathlib import Path
from typing import Union
from tqdm import tqdm

class ImportLammpsDumpposes():
    def __init__(self):
        pass

    def import_lammps_dumppos(self, ifn: str):
        '''
        LAMMPSで作成されたdumpposファイルを読み込む(1つのdumppos)
        '''

        with open(ifn, 'r') as ifp:
            # "ITEM: TIMESTEP"
            spline = ifp.readline().split()
            time_step = int(ifp.readline())
            
            # ITEM: NUMBER OF ATOMS
            spline = ifp.readline().split()
            total_atoms = int(ifp.readline())

            # ITEM: BOX BOUNDS xy xz yz pp pp pp
            spline = ifp.readline().split()
            for dim in range(3):
                spline = ifp.readline().split()
                self.cell[dim] = float(spline[1])
            
            # ITEM: ATOMS id type x y z
            spline = ifp.readline().split()
            columns = spline[2:]

        max_column_num = max(9, len(columns) + 1)
        df = pd.read_csv(ifn, sep='\s+', names=range(max_column_num))

        self.step_nums = [None for _ in range(len(df) // (total_atoms + 9))]
        self.step_num_to_step_idx = {}
        self.atoms = [None for _ in range(len(self.step_nums))] 

        for step_idx in range(len(self.step_nums)):
            self.atoms[step_idx] = df.iloc[step_idx * (total_atoms + 9):(step_idx + 1) * (total_atoms + 9),:].copy()
            step_num = self.atoms[step_idx].iloc[1, 0]
            step_num = int(step_num)
            self.step_nums[step_idx] = step_num
            self.step_num_to_step_idx[step_num] = step_idx
            self.atoms[step_idx] = self.atoms[step_idx].iloc[9:,:len(columns)]
            self.atoms[step_idx].columns = columns
            self.atoms[step_idx]['id'] = self.atoms[step_idx]['id'].astype(int)
            self.atoms[step_idx]['type'] = self.atoms[step_idx]['type'].astype(int)
            self.atoms[step_idx][['x', 'y', 'z']] = self.atoms[step_idx][['x', 'y', 'z']].astype(float)
            self.atoms[step_idx].sort_values('id', inplace=True)
            # 0-index
            self.atoms[step_idx].index = self.atoms[step_idx]['id'] - 1
            self.atoms[step_idx].drop('id', axis=1, inplace=True)

    def import_lammps_dumpposes(self, input_file_dir: Union[str, Path]='./', skip_num=1):
        '''
        LAMMPSで作成されたdumpposファイルを読み込む(複数のdumppos)
        '''
        input_file_dir = Path(input_file_dir)
        dumppos_file_paths = list(input_file_dir.glob(f'./dump.pos.*'))
        step_nums = []
        for dumppos_file_path in dumppos_file_paths:
            step_nums.append(int(dumppos_file_path.name[9:]))
        step_nums.sort()
        step_nums = step_nums[::skip_num]
        self.step_nums = step_nums
        self.step_num_to_step_idx = dict()
        self.atoms = [None] * len(step_nums)
        
        for step_idx, step_num in enumerate(tqdm(step_nums, desc='[importing lammps dumpposes]')):
            dumppos_file_path = f'dump.pos.{step_num}'
            with open(dumppos_file_path, 'r') as ifp:
                # "ITEM: TIMESTEP"
                spline = ifp.readline().split()
                step_num = int(ifp.readline())
                self.step_nums[step_idx] = step_num
                self.step_num_to_step_idx[step_num] = step_idx
                # ITEM: NUMBER OF ATOMS
                spline = ifp.readline().split()
                total_atoms = int(ifp.readline())

                # ITEM: BOX BOUNDS xy xz yz pp pp pp
                spline = ifp.readline().split()
                for dim in range(3):
                    spline = ifp.readline().split()
                    self.cell[dim] = float(spline[1])
                
                # ITEM: ATOMS id type x y z
                spline = ifp.readline().split()
                columns = spline[2:]

            df = pd.read_csv(dumppos_file_path, sep='\s+', names=columns, skiprows=9)
            df[['x', 'y', 'z']] = df[['x', 'y', 'z']].astype(float)
            df['type'] = df['type'].astype(int)
            df['id'] = df['id'].astype(int)
            df.sort_values('id', inplace=True)
            df.index = df['id'] - 1
            df.drop('id', axis=1, inplace=True)

            if 'vx' in df.columns and 'vy' in df.columns and 'vz' in df.columns:
                df[['vx', 'vy', 'vz']] = df[['vx', 'vy', 'vz']].astype(float)
            
            self.atoms[step_idx] = df
