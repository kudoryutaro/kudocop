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


class AnalyzeBondForSDats():
    def __init__():
        pass

    def count_bonds(self, cut_off) -> pd.DataFrame:
        self.get_connect_lists(cut_off)

        bond_types = []
        for atom_type1 in self.atom_type_set:
            for atom_type2 in self.atom_type_set:
                if (atom_type1, atom_type2) in bond_types or (atom_type2, atom_type1) in bond_types:
                    continue
                if atom_type2 < atom_type2:
                    continue
                bond_types.append((atom_type1, atom_type2))

        bond_counter = [{bond_type: 0 for bond_type in bond_types}
                        for _ in range(len(self.step_nums))]

        for step_idx, step_num in enumerate(trange(len(self.step_nums), desc='[counting bonds]')):

            self_atom_type = self.atoms[step_idx]['type'].values
            for atom_idx, c_list in enumerate(self.connect_lists[step_idx]):
                atom_type = self_atom_type[atom_idx]
                for next_atom_idx in c_list:
                    next_atom_type = self_atom_type[next_atom_idx]
                    if next_atom_type < atom_type:
                        continue
                    bond_counter[step_idx][(atom_type, next_atom_type)] += 1

            for (atom_type, next_atom_type) in bond_counter[step_idx]:
                if atom_type == next_atom_type:
                    bond_counter[step_idx][(atom_type, next_atom_type)] //= 2
        df_bond_count = pd.DataFrame(bond_counter, index=self.step_nums)
        print(len(df_bond_count))

        def change_column_name(column):
            atom_type1, atom_type2 = column
            return f'{self.atom_type_to_symbol[atom_type1]}-{self.atom_type_to_symbol[atom_type2]}'

        df_bond_count = df_bond_count.rename(columns=change_column_name)

        return df_bond_count
