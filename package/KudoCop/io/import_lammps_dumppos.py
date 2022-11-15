import numpy as np
import pandas as pd
import sys
from pathlib import Path

class ImportLammpsDumppos():
    def __init__(self):
        pass

    def import_lammps_dumppos(self, ifn: str):
        '''
        LAMMPSで作成されたdumpposファイルを読み込む
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

        df = pd.read_csv(ifn, sep='\s+', names=columns, skiprows=9)
        df[['x', 'y', 'z']] = df[['x', 'y', 'z']].astype(float)
        df['type'] = df['type'].astype(int)
        df['id'] = df['id'].astype(int)
        df.sort_values('id', inplace=True)
        df.index = df['id'] - 1
        df.drop('id', axis=1, inplace=True)

        if 'vx' in df.columns and 'vy' in df.columns and 'vz' in df.columns:
            df[['vx', 'vy', 'vz']] = df[['vx', 'vy', 'vz']].astype(float)
        
        self.atoms = df
