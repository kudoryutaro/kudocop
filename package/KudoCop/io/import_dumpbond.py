import numpy as np
import pandas as pd


class ImportDumpbond():
    def __init__(self):
        pass

    def import_dumpbond(self, ifn: str) -> None:
        with open(ifn, 'r') as ifp:
            current_row = 0
            while True:
                current_row += 1
                spline = ifp.readline().split()
                if len(spline) == 0:
                    continue
                if spline[0] == 'ITEM:' and spline[1] == 'NUMBER':
                    total_atoms = int(ifp.readline())
                    current_row += 1
                    break
        skip_rows = current_row
        df_bond = pd.read_csv(ifn, skiprows=skip_rows,
                              sep=' ', names=[0, 1, 2, 3, 4])
        df_bond = df_bond.replace('Atom', -1)
        df_bond[0] = df_bond[0].astype(int)
        atom_ids_and_neibour_nums = df_bond.iloc[:, [
            1, 2]][df_bond.iloc[:, 3].isna()].astype(int)

        self.bondorder_list = [None] * total_atoms
        self.bondorder_connect_list = [None] * total_atoms

        df_bond_order_values = df_bond.iloc[:, [1, 2, 3, 4]].values

        # 0-indexed
        df_bond_neibour_values = df_bond.iloc[:, 0].values - 1
        for idx, atom_id, neibour_num in atom_ids_and_neibour_nums.itertuples():

            # atom_id - 1 == atom_idx
            self.bondorder_list[atom_id -
                                1] = df_bond_order_values[idx + 1:idx + 1 + neibour_num]
            self.bondorder_connect_list[atom_id -
                                        1] = df_bond_neibour_values[idx + 1:idx + 1 + neibour_num]
