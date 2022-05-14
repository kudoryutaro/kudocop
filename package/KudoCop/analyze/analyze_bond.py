import sys
from tqdm import tqdm, trange
import pandas as pd
import numpy as np

class AnalyzeBond():
    def __init__():
        pass
    

    def count_bonds(self, cut_off) -> dict:
        self.get_connect_list(cut_off)

        # bond_counter[(atom_type1,atom_type2)] atom1,atom2の結合数
        bond_counter = dict()

        for atom_type1 in self.atom_type_set:
            for atom_type2 in self.atom_type_set:
                if (atom_type1, atom_type2) in bond_counter or (atom_type2, atom_type1) in bond_counter:
                    continue
                if atom_type2 < atom_type2:
                    continue
                bond_counter[(atom_type1, atom_type2)] = 0

        self_atom_type = self.atoms['type'].values
        for atom_idx, c_list in enumerate(self.connect_list):
            atom_type = self_atom_type[atom_idx]
            for next_atom_idx in c_list:
                next_atom_type = self_atom_type[next_atom_idx]
                if next_atom_type < atom_type:
                    continue
                bond_counter[(atom_type, next_atom_type)] += 1

        for (atom_type, next_atom_type) in bond_counter:
            if atom_type == next_atom_type:
                bond_counter[(atom_type, next_atom_type)] //= 2

        return bond_counter


    # def count_bonds_selfs(selfs, cut_off):
    #     dfs_count_bonds = []
    #     for step_idx, step_num in enumerate(tqdm(selfs.step_nums)):
    #         counter = selfs.data[step_idx].count_bonds(cut_off)
    #         dfs_count_bonds.append(pd.DataFrame(counter, [step_num]))
    #     df_count_bonds = pd.concat(dfs_count_bonds).fillna(0)
    #     df_count_bonds.columns = list(df_count_bonds.columns)

    #     def change_column_name(column):
    #         atom_type1, atom_type2 = column
    #         return f'{selfs.data[0].atom_type_to_symbol[atom_type1]}-{selfs.data[0].atom_type_to_symbol[atom_type2]}'

    #     df_count_bonds = df_count_bonds.rename(columns=change_column_name)

    #     return df_count_bonds
