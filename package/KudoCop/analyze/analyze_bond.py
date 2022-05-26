import sys
from tqdm import tqdm, trange
import pandas as pd
import numpy as np
import itertools
import math

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
                if atom_type2 < atom_type1:
                    continue
                bond_counter[(atom_type1, atom_type2)] = 0

        sdat_atom_type = self.atoms['type'].values
        for atom_idx, c_list in enumerate(self.connect_list):
            atom_type = sdat_atom_type[atom_idx]
            for next_atom_idx in c_list:
                next_atom_type = sdat_atom_type[next_atom_idx]
                if next_atom_type < atom_type:
                    continue
                bond_counter[(atom_type, next_atom_type)] += 1

        for (atom_type, next_atom_type) in bond_counter:
            if atom_type == next_atom_type:
                bond_counter[(atom_type, next_atom_type)] //= 2

        return bond_counter

    def get_bond_angle(self, cut_off) -> list:
        self.get_connect_list(cut_off)
        max_atom_type = max(self.atom_type_set)
        # bond_angle[neibour1][mid][neibour2] : neibour1-mid-neibour2 の角度の入ったリスト
        bond_angle = [[[[] for _ in range(max_atom_type + 1)] for _ in range(
            max_atom_type + 1)] for _ in range(max_atom_type + 1)]
        sdat_atom_type = self.atoms['type'].values
        sdat_atom_xyz = self.atoms[['x', 'y', 'z']].values
        for mid_atom_idx, c_list in enumerate(self.connect_list):
            mid_atom_type = sdat_atom_type[mid_atom_idx]
            for neibour_atom_idx1, neibour_atom_idx2 in itertools.combinations(c_list, 2):
                neibour_atom_type1 = sdat_atom_type[neibour_atom_idx1]
                neibour_atom_type2 = sdat_atom_type[neibour_atom_idx2]
                angle = calc_angle_of_ABC(
                    sdat_atom_xyz[neibour_atom_idx1], sdat_atom_xyz[mid_atom_idx], sdat_atom_xyz[neibour_atom_idx2])
                if neibour_atom_type1 == neibour_atom_type2:
                    bond_angle[neibour_atom_type1][mid_atom_type][neibour_atom_type2].append(
                        angle)
                else:
                    bond_angle[neibour_atom_type1][mid_atom_type][neibour_atom_type2].append(
                        angle)
                    bond_angle[neibour_atom_type2][mid_atom_type][neibour_atom_type1].append(
                        angle)
        return bond_angle

    def get_coordination_number(self, cut_off):
        self.get_connect_list(cut_off)
        bond_counter = dict()

        for atom_type1 in self.atom_type_set:
            for atom_type2 in self.atom_type_set:
                bond_counter[(atom_type1, atom_type2)] = 0

        sdat_atom_type = self.atoms['type'].values
        for atom_idx, c_list in enumerate(self.connect_list):
            atom_type = sdat_atom_type[atom_idx]
            for next_atom_idx in c_list:
                next_atom_type = sdat_atom_type[next_atom_idx]
                bond_counter[(atom_type, next_atom_type)] += 1

        for (atom_type, next_atom_type) in bond_counter:
            if atom_type == next_atom_type:
                bond_counter[(atom_type, next_atom_type)] //= 2
            atom_number = len(self.atoms[self.atoms['type'] == atom_type])
            if atom_number != 0:
                bond_counter[(atom_type, next_atom_type)] /= atom_number

        return bond_counter

    # def get_bond_length(self, cut_off):
    #     self.get_connect_list(cut_off)
    #     bond_length = dict()

    #     for atom_type1 in self.atom_type_set:
    #         for atom_type2 in self.atom_type_set:
    #             if (atom_type1, atom_type2) in bond_length or (atom_type2, atom_type1) in bond_length:
    #                 continue
    #             if atom_type2 < atom_type1:
    #                 continue
    #             bond_length[(atom_type1, atom_type2)] = []

    #     sdat_atom_type = self.atoms['type'].values
    #     sdat_atom_xyz = self.atoms[['x', 'y', 'z']].values


    #     for atom_idx, c_list in enumerate(self.connect_list):
    #         atom_type = sdat_atom_type[atom_idx]
    #         atom_xyz = sdat_atom_xyz[atom_idx]
    #         for next_atom_idx in c_list:
    #             next_atom_type = sdat_atom_type[next_atom_idx]
    #             next_atom_xyz = sdat_atom_xyz[next_atom_idx]
    #             if next_atom_type < atom_type:
    #                 continue
    #             bond_length[(atom_type, next_atom_type)
    #                             ].append(math.sqrt(sum((atom_xyz - next_atom_xyz) ** 2)))
    #     bond_counter = self.count_bonds(cut_off)
    #     for (atom_type, next_atom_type) in bond_length:
    #         if atom_type == next_atom_type:
    #             bond_length[(atom_type, next_atom_type)] /= 2
    #         if bond_counter[(atom_type, next_atom_type)] != 0:
    #             bond_length[(atom_type, next_atom_type)
    #                             ] /= bond_counter[(atom_type, next_atom_type)]

    #     return bond_length
    

        



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
                if atom_type2 < atom_type1:
                    continue
                bond_types.append((atom_type1, atom_type2))

        bond_counter = [{bond_type: 0 for bond_type in bond_types}
                        for _ in range(len(self.step_nums))]

        for step_idx, step_num in enumerate(trange(len(self.step_nums), desc='[counting bonds]')):

            sdat_atom_type = self.atoms[step_idx]['type'].values
            for atom_idx, c_list in enumerate(self.connect_lists[step_idx]):
                atom_type = sdat_atom_type[atom_idx]
                for next_atom_idx in c_list:
                    next_atom_type = sdat_atom_type[next_atom_idx]
                    if next_atom_type < atom_type:
                        continue
                    bond_counter[step_idx][(atom_type, next_atom_type)] += 1

            for (atom_type, next_atom_type) in bond_counter[step_idx]:
                if atom_type == next_atom_type:
                    bond_counter[step_idx][(atom_type, next_atom_type)] //= 2
        df_bond_count = pd.DataFrame(bond_counter, index=self.step_nums)

        def change_column_name(column):
            atom_type1, atom_type2 = column
            return f'{self.atom_type_to_symbol[atom_type1]}-{self.atom_type_to_symbol[atom_type2]}'

        df_bond_count = df_bond_count.rename(columns=change_column_name)

        return df_bond_count


def calc_angle_of_ABC(a: np.array, b: np.array, c: np.array) -> float:
    """
    ABとBCのなす角度(degree)を返す
    """
    # ベクトルを定義
    vec_a = a - b
    vec_c = c - b

    # コサインの計算
    length_vec_a = np.linalg.norm(vec_a)
    length_vec_c = np.linalg.norm(vec_c)
    inner_product = np.inner(vec_a, vec_c)
    cos = inner_product / (length_vec_a * length_vec_c)

    # 角度（ラジアン）の計算
    rad = np.arccos(cos)

    # 弧度法から度数法（rad ➔ 度）への変換
    degree = np.rad2deg(rad)

    return degree
