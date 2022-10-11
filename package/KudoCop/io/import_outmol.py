import numpy as np
import pandas as pd
import sys
from pathlib import Path

class ImportOutmol():
    def __init__(self):
        pass

    def import_outmol_opt(self, ifn: str):
        '''
        構造最適化によって最適化された最後の構造を読み込む
        '''
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for index, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if len(spline) >= 2 and spline[0] == 'Final' and spline[1] == 'Coordinates':
                start_index = index + 3
            if len(spline) >= 2 and spline[0] == '$cell' and spline[1] == 'vectors':
                self.cell[0] = round(float(lines[index + 1].split()
                                     [0]) / 1.889726124993590, 2)
                self.cell[1] = round(float(lines[index + 2].split()
                                     [1]) / 1.889726124993590, 2)
                self.cell[2] = round(float(lines[index + 3].split()
                                     [2]) / 1.889726124993590, 2)

        index = start_index
        data = []

        while True:
            spline = lines[index].split()
            try:
                int(spline[0])
                data.append(spline)
            except:
                break
            index += 1

        df_outmol = pd.DataFrame(data, columns=['id', 'type', 'x', 'y', 'z'])
        df_outmol.drop('id', axis=1, inplace=True)
        df_outmol[['x', 'y', 'z']] = df_outmol[['x', 'y', 'z']].astype(float)

        df_outmol['type'] = df_outmol['type'].replace(self.atom_symbol_to_type)

        self.atoms = df_outmol.copy()


class ImportOutmolForSDats():
    def __init__():
        pass

    def import_outmol_md(self, ifn: Path):
        """dmol3で第一原理分子動力学をした結果のoutmolファイルを読み込む
        Parameters
        ----------
            ifn : str or Path
                outmolファイルのパス
        Note
        ----
            step 0 は読み込まない
            step (0, number_of_step]を読み込む
            (stepの数をnumber_of_step個にしたいため)
        """
        with open(ifn, 'r') as f:
            lines = f.readlines()
            splines = list(map(lambda l:l.split(), lines))
        
        # read cell size , total_atoms, number_of_steps
        for spline_idx, spline in enumerate(splines):
            if len(spline) == 0:
                continue
            if spline[0] == '$cell':
                self.cell = [None, None, None]
                for i in range(3):
                    self.cell[i] = float(splines[spline_idx + 1 + i][i])
            if spline[0] == 'N_atoms':
                total_atoms = int(spline[2])
            if spline[0] == 'Step':
                number_of_steps = int(splines[spline_idx + 1][1])
        
        self.step_nums = list(range(1, number_of_steps + 1))
        self.step_num_to_step_idx = {
            step_num: step_idx for step_idx, step_num in enumerate(self.step_nums)
        }

        self.atoms = [None] * len(self.step_nums)
        self.force = [None] * len(self.step_nums)
        self.potential_energy = [None] * len(self.step_nums) # Unit : eV

        # read atoms
        current_step_idx = 0
        current_step_num = 0
        for spline_idx, spline in enumerate(splines):
            if current_step_idx > len(self.step_nums):
                break
            if len(spline) >= 3 and spline[0] == 'df' and spline[1] == 'ATOMIC' and spline[2] == 'COORDINATES':
                if current_step_num == 0:
                    current_step_num += 1
                    continue
                atoms_line = []
                for i in range(2, total_atoms + 2):
                    spline_atom_idx = spline_idx + i
                    atoms_line.append(splines[spline_atom_idx][1:5])
                df_atoms = pd.DataFrame(data=atoms_line, columns=['type', 'x', 'y', 'z'])
                df_atoms['type'] = df_atoms['type'].map(self.atom_symbol_to_type)
                self.atoms[current_step_idx] = df_atoms
                current_step_idx += 1
                current_step_num += 1
        
        # read velocity and acceleration
        current_step_idx = 0
        current_step_num = 0
        for spline_idx, spline in enumerate(splines):
            if current_step_idx > len(self.step_nums):
                break
            if len(spline) >= 3 and spline[0] == 'dq' and spline[1] == 'ATOMIC' and spline[2] == 'VELOCITIES':
                if current_step_num == 0:
                    current_step_num += 1
                    continue
                atoms_line = []
                for i in range(2, total_atoms + 2):
                    spline_atom_idx = spline_idx + i
                    atoms_line.append(splines[spline_atom_idx][2:8])
                self.atoms[current_step_idx][['vx', 'vy', 'vz', 'ax', 'ay', 'az']] = atoms_line
                current_step_idx += 1
                current_step_num += 1

        # read potential energy
        current_step_idx = 0
        for spline_idx, spline in enumerate(splines):
            if current_step_idx >= len(self.step_nums):
                break
            # if len(spline) >= 5 and spline[0] == 'Step' and spline[1] == 'Kin.+Pot.' and spline[2] == 'Energy' and spline[3] == 'Pot.' and spline[4] == 'Energy':
            #     self.potential_energy[current_step_idx] = float(splines[spline_idx + 1][4]) * 27.2114
            #     current_step_idx += 1
            if len(spline) >= 5 and spline[0] == 'Step' and spline[1] == 'System' and spline[2] == 'Energy' and \
                spline[3] == 'Pot.' and spline[4] == 'Energy':
                self.potential_energy[current_step_idx] = float(splines[spline_idx + 1][4]) * 27.2114
                current_step_idx += 1

