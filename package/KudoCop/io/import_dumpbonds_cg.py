import numpy as np
import pandas as pd
from tqdm import tqdm

class ImportDumpbondsCG():
    def __init__(self):
        pass

    def import_dumpbonds_cg(self) -> None:
        
        for step_idx, step_num in enumerate(tqdm(self.step_nums, desc='[importing dumpbond]')):
            current_file_name = f'{self.dir_name}/dump.bond.{step_num}'
            with open(current_file_name, 'r') as ifp:
                lines = ifp.readlines()
            total_atoms =len(lines) -1
            self.connect_list = [[] for _ in range(total_atoms)]
            for line in lines[1:]:
                spline = list(map(int,line.split()))
                atom_idx = spline[0] - 1
                c_list = spline[2:]
                for i in range(len(c_list)):
                    c_list[i] -= 1

                self.connect_list[atom_idx] = c_list


            for step_idx, step_num in enumerate(tqdm(self.step_nums, desc='[importing dumpbond]')):
                current_file_name = f'{self.dir_name}/dump.bond.{step_num}'
                df_bond = pd.read_csv(current_file_name, skiprows=skip_rows,
                                    sep=' ', names=[0, 1, 2, 3, 4])
                df_bond = df_bond.replace('Atom', -1)
                df_bond[0] = df_bond[0].astype(int)
                atom_ids_and_neibour_nums = df_bond.iloc[:, [
                    1, 2]][df_bond.iloc[:, 3].isna()].astype(int)

                df_bond_order_values = df_bond.iloc[:, [1, 2, 3, 4]].values

                # 0-indexed
                df_bond_neibour_values = df_bond.iloc[:, 0].values - 1
                for idx, atom_id, neibour_num in atom_ids_and_neibour_nums.itertuples():

                    # atom_id - 1 == atom_idx
                    self.bondorder_lists[step_idx][atom_id -
                                                1] = df_bond_order_values[idx + 1:idx + 1 + neibour_num]
                    self.bondorder_connect_lists[step_idx][atom_id -
                                                        1] = df_bond_neibour_values[idx + 1:idx + 1 + neibour_num]
