import numpy as np
import pandas as pd
from tqdm import tqdm


class ImportDumpbondsCG():
    def __init__(self):
        pass

    def import_dumpbonds_cg(self, max_coordination_num) -> None:

        self.connect_lists_from_dumpbond_cg = [
            None for _ in range(len(self.step_nums))]
        for step_idx, step_num in enumerate(tqdm(self.step_nums, desc='[importing dumpbond]')):
            current_file_name = f'{self.dir_name}/dump.bond.{step_num}'

            df = pd.read_csv(current_file_name, skiprows=1, sep=' ',
                             header=None, names=range(max_coordination_num))
            self.connect_lists_from_dumpbond_cg[step_idx] = [
                None for _ in range(len(df.index))]
            df.drop([0, 2], axis=1, inplace=True)
            df = df.fillna(-1).astype(int) - 1

            values = df.values

            for row in values:
                atom_idx = row[0]
                for col_i, next_atom_idx in enumerate(row):
                    if next_atom_idx < 0:
                        break
                self.connect_lists_from_dumpbond_cg[step_idx][atom_idx] = row[1:col_i]
