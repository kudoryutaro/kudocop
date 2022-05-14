import numpy as np
import pandas as pd
from tqdm import tqdm


class ImportDumpposes():
    def __init__(self):
        pass

    def import_dumpposes(self) -> None:
        current_row = 0
        first_step_file_name = f'{self.dir_name}/dump.pos.{self.step_nums[0]}'
        with open(first_step_file_name, 'r') as ifp:
            while True:
                current_row += 1
                spline = ifp.readline().split()
                if len(spline) == 0:
                    continue
                if spline[0] == "ITEM:" and spline[1] == "BOX":
                    for dim in range(3):
                        spline = ifp.readline().split()
                        self.cell[dim] = float(spline[1])
                    current_row += 3
                    continue
                if spline[0] == 'ITEM:' and spline[1] == 'ATOMS':
                    columns = spline[3:]
                    break

        skip_rows = current_row

        for step_idx, step_num in enumerate(tqdm(self.step_nums, desc='[importing dumppos]')):
            current_file_name = f'{self.dir_name}/dump.pos.{step_num}'
            self.atoms[step_idx] = pd.read_csv(
                current_file_name, skiprows=skip_rows, sep=' ', names=columns)
            self.atoms[step_idx][['type', 'mask']
                                 ] = self.atoms[step_idx][['type', 'mask']].astype(int)
            self.atoms[step_idx].index = self.atoms[step_idx].index - 1
            self.atoms[step_idx].sort_index(inplace=True)
            self.atom_type_set |= set(self.atoms[step_idx]['type'])
